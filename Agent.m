%% A class for defining the agents (robots) in the MRTA class
% 
% The class contains the properties of the agent, such as the tasks it can
% perform, as well as methods that are useful for building the MRTA problem.
% 

classdef Agent < handle

    properties
        Gv % Location graph of the agent
        name % Name of the agent
        loc % Location of the agent
        l_with_t % Array containing the indices of locations where the agent can perform at least one task
        n_l_with_t % Number of locations where the agent can perform at least one task
        ts_at_l % Array containing the number of tasks the agent can perform at each location
        t_loc % Local copy of the current pending tasks
        Eat % Array defining the tasks it can perform (row of the global Eat for this agent)
        Evt % Matrix defining the tasks it can perform in each location (takes into account if the task can be assigned to that location)
        idx_Evt % Stores the indeces of each location-task pair in the dv part used to enconde those tasks
        task_idx % Stores the indeces of the locations and tasks that the agent can do (the rows and cols of Evt that are 1)
        dv % Struct defining the dimensions associated with its decision variables
           % dv.idx - Integer that identifies the index in the decision varibales corresponding to this agent
           % dv.n_e - Number of edges that the agent can travel through
           % dv.n - Number of locations where the agent can stop
           % dv.n_t - Number of tasks-location the agent can do
           % dv.nz - Total number of decision variables associated with the agent

        % For defining and working with its own MILP problem
        miqp % Instance of MIQPholder that defines the agent's MILP problem
        % NOTE: I state the problem as an MIQP for generality, although we use an MILP to pose the MRTA problem
        param % Structure containing information like the cost of each task, etc.

    end
    
    methods

        % Constructor
        function obj = Agent(Gv, Evt, Eat, opt)
            arguments
                Gv digraph
                Evt (:, :) % Matrix defining the tasks it can perform in each location (takes into account if the task can be assigned to that location)
                Eat (:, :) % Array defining the tasks it can perform (row of the global Eat for this agent)
                opt.name string = "" % Name of the agent
                opt.loc (1, 1) = 0 % Current location of the agent
                opt.t_loc = {} % Local copy of the current pending tasks
                opt.param struct = struct() % Contains parameters of the MRTA problem, such as costs
                opt.idx (1, 1) = 0 % Index where the agent starts in the global decision vector
            end
            obj.Gv = Gv;
            obj.name = opt.name;
            obj.param = opt.param;
            obj.loc = opt.loc;
            obj.t_loc = opt.t_loc;
            obj.Eat = Eat;
            obj.Evt = Eat.*Evt;
            obj.dv.idx = opt.idx;
            obj.dv.n_e = size(obj.Gv.Edges, 1);
            obj.dv.n = size(obj.Gv.Nodes, 1);
            obj.dv.n_t = sum(Eat*Evt');
            obj.dv.nz = obj.dv.n_e + 2*obj.dv.n + obj.dv.n_t;
            
            curr_idx_Evt = 1;
            obj.idx_Evt = obj.Evt;
            dim_Evt = size(Evt);
            for i = 1:dim_Evt(1)
                for j = 1:dim_Evt(2)
                    if obj.idx_Evt(i, j) ~= 0
                        obj.idx_Evt(i, j) = curr_idx_Evt;
                        curr_idx_Evt = curr_idx_Evt + 1;
                    end
                end
            end

            obj.ts_at_l = sum(obj.Evt, 2);
            obj.n_l_with_t = nnz(obj.ts_at_l);
            obj.l_with_t = zeros(1, obj.n_l_with_t);
            j = 1;
            for i = 1:size(obj.Evt, 1)
                if obj.ts_at_l(i) > 0
                    obj.l_with_t(j) = i;
                    j = j+1;
                end
            end

            [obj.task_idx.v, obj.task_idx.t] = find(obj.Evt == 1);

        end

        function build_MILP(self)
            % Computes the MILP problem of the agent (without taking into account single-task constraints between agents)
            % Stores the MILP problem in self.miqp
            
            % Cost function
            c = self.build_c();

            % Equality constraints
            [Cio, dio] = self.gen_in_out_eq_constraints();
            % [Css, dss] = self.gen_single_stop_eq_constraints();

            % Inequality constraints
            [Ccd, lcd, ucd] = self.gen_can_do_ineq_constraints();
            [Cst, lst, ust] = self.gen_subtour_ineq_constraints();

            % Box constraints
            [lbz, ubz] = self.gen_box_constraints();

            % Integral constraints
            [D, li, ui] = self.gen_integer_constraints();

            % Construct MILP problem
            Q = zeros(self.dv.nz); % We consider a MILP problem in this case
            A = [Ccd; Cst];
            lb = [lcd; lst];
            ub = [ucd; ust];
            C = Cio;
            b = dio;
            idx_int = self.get_idx_integrality();
            self.miqp = MIQPholder(Q, c, A, lb, ub, C, b, D, li, ui, lbz, ubz, idx_int);

        end
        
        function c = build_c(self)
            % Returns the basic c vector with the costs associated to the decision variables of this agent
            c = zeros(self.dv.nz, 1);
            c(1:self.dv.n_e) = self.Gv.Edges.Weight; % Add the edge weights (move between locations)
            c(self.dv.n_e + (1:self.dv.n)) = self.param.cost_stop; % Add stopping cost in a location
            w_tasks = (self.param.cost_task.*self.Evt)'; % Add cost of completing task
            w_tasks = w_tasks(self.Evt' == 1); % Extract elements in the dv vector 
            c(self.dv.n_e + self.dv.n + (1:self.dv.n_t)) = w_tasks; % Add these costs to c
            c(self.dv.n_e + self.dv.n + self.dv.n_t + (1:self.dv.n)) = self.param.cost_times; % Add cost for the time-variables
        end

        function [Cio, dio] = gen_in_out_eq_constraints(self)
            % This function returns matrix Cio and vector dio for the equality
            % constraints that force in=out at each location (considering starting location)
            Cio = zeros(self.dv.n, self.dv.nz);
            for l = 1:self.dv.n
                Cio(l, self.Gv.inedges(l)) = -1;
                Cio(l, self.Gv.outedges(l)) = 1;
                Cio(l, self.dv.n_e + l) = 1;
            end
            dio = self.update_dio();
        end

        function dio = update_dio(self, in_loc)
            % Updates the vector dio of gen_in_out_eq_constraints for the current agent location @in_loc.
            % If no @in_loc is provided, then the agent't self.loc variable is used.
            % If an @in_loc is provided, then the agent updates its self.loc to the new location.          
            arguments
                self Agent
                in_loc (1, 1) = NaN
            end
            if isnan(in_loc)
                in_loc = self.loc;
            else
                self.loc = in_loc;
            end
            dio = zeros(self.dv.n, 1);
            dio(in_loc) = 1;
        end

        function [Css, dss] = gen_single_stop_eq_constraints(self)
            % This function returns matrix Css and vector dxx for the equality
            % constraints that force an agent to only stop once
            Css = zeros(1, self.dv.nz);
            Css(1, self.dv.n_e + (1:self.dv.n)) = ones(1, self.dv.n);
            dss = 1;
        end

        function [Ccd, lcd, ucd] = gen_can_do_ineq_constraints(self)
            % This function returns matrix Ccd and vectos lcd, ucd for the inequality
            % constraints that determine if an agent can do a task (only if it has visited the location): lcd <= Ccd*z <= ucd
            Ccd = zeros(self.n_l_with_t, self.dv.nz);
            v = 1;
            k = size(self.Evt, 2);
            for l = self.l_with_t
                Ccd(v, self.Gv.inedges(l)) = -1;
                for t = 1:k
                    [has_task, ~, idx_tv] = self.has_task_loc(t, l);
                    if has_task
                        Ccd(v, idx_tv) = 1/self.ts_at_l(l);
                    end
                end
                v = v+1;
            end
            [lcd, ucd] = self.update_cd_bounds();
        end

        function [lcd, ucd] = update_cd_bounds(self, in_loc)
            % Updates the lower and upper bounds lcd and ucd of gen_can_do_ineq_constraints() for the current agent location @in_loc.
            % If no @in_loc is provided, then the agent's self.loc variable is used.
            % If an @in_loc is provided, then the agent updates its self.loc to the new location.          
            arguments
                self Agent
                in_loc (1, 1) = NaN
            end
            if isnan(in_loc)
                in_loc = self.loc;
            else
                self.loc = in_loc;
            end
            ucd = zeros(self.n_l_with_t, 1);
            if self.ts_at_l(in_loc) > 0
                ucd(self.l_with_t == in_loc) = 1;
            end
            lcd = -inf*ones(self.n_l_with_t, 1);
        end

        function [Cst, lst, ust] = gen_subtour_ineq_constraints(self, opt)
            % This function returns matrix Cst and vectos lst, ust for the inequality
            % constraints that avoid subtours in the optimal solution: lst <= Cst*z <= ust
            arguments
                self Agent
                opt.with_box (1, 1) {mustBeNumericOrLogical} = false;
            end

            W = sum(self.Gv.Edges.Weight); % Total travel distance in the location graph. This is a safe upper bound for the soubtour variables
            
            % Constraints to impose monotonic increase whenever the agent goes from one location to another
            lst = -W*ones(self.dv.n_e, 1);
            ust = zeros(self.dv.n_e, 1);

            Cst = zeros(self.dv.n_e, self.dv.nz);
            idx_t_var = self.dv.n_e + self.dv.n + self.dv.n_t;
            for edge = 1:self.dv.n_e
                [i, j] = findedge(self.Gv, edge);
                Cst(edge, edge) = W;
                Cst(edge, idx_t_var+i) = 1;
                Cst(edge, idx_t_var+j) = -1;
                ust(edge) = W - self.Gv.Edges.Weight(edge);
            end

            if opt.with_box
                Cst = [Cst; zeros(self.dv.n, self.dv.nz-self.dv.n), eye(self.dv.n)];
                lst = [lst; zeros(self.dv.n, 1)];
                ust = [ust; W*ones(self.dv.n, 1)];
            end

        end

        function [lbz, ubz] = gen_box_constraints(self)
            % This function returns the box constraints lbz and ubz for z: lbz <= z <= ubz
            lbz = zeros(self.dv.nz, 1);
            ubz = [ones(self.dv.n_e + self.dv.n + self.dv.n_t, 1); self.dv.n*ones(self.dv.n, 1)];
            % NOTE: the last part of ubz is for the subtour variables. The bound might change if a different step is taken between them
            %       For now dv.n is enough, because times are always normalized to 1, so accumulated time will never be larger than dv.n.
        end

        function [D, li, ui] = gen_integer_constraints(self)
            % This function returns matrix D and vectos li, ui that determine the integrality constraints
            D = [eye(self.dv.n_e + self.dv.n + self.dv.n_t), zeros(self.dv.n_e + self.dv.n + self.dv.n_t, self.dv.n)];
            li = zeros(self.dv.n_e + self.dv.n + self.dv.n_t, 1);
            ui = ones(self.dv.n_e + self.dv.n + self.dv.n_t, 1);
        end

        function int_idx = get_idx_integrality(self)
            % This function returns a vector containing the indices of the MILP problem that are integer variables
            int_idx = 1:(self.dv.n_e + self.dv.n + self.dv.n_t); % NOTE: Assuming here that the subtour varibales are the last in z
        end

        function [has, idx, idx_local] = has_task_loc(self, t, v)
            % Returns true if the agent has the given task @t at location @v
            % Returns the index of the decision variables associated to the task-location pair.
            % Also returns the local index (instead of the one related to the complete dv vector).
            % Returns an empty [] if the agent does not have that task associated to that location
            if ~( self.idx_Evt(v, t) == 0 )
                has = true;
                idx_local = self.dv.n_e + self.dv.n + self.idx_Evt(v, t);
                idx = idx_local + self.dv.idx;
            else
                has = false;
                idx_local = [];
                idx = [];
            end
        end

        function z = round_z(self, z)
            % Returns a vector @z_round where the integral elements of @z have been rounded to 0 or 1
            for i = 1:self.dv.n_e + self.dv.n + self.dv.n_t
                if z(i) >= 0.5
                    z(i) = 1;
                else
                    z(i) = 0;
                end
            end

        end

        function text = res2tex(self, z, Gv, names)
            % Generates the text that iterprets the solution @z corresponding to this agent.
            % To print human-interpretable text, the function requires the location graph @Gv
            % and the @names structure containing the names of locations and tasks
            % @z is assumed to be rounded so that it contains numbers with up to one significant digit (and integers 1 or 0 are exact integers)

            % Print basic information
            text = "";
            text = strcat(text, "Agent ", self.name, "\n");
            text = strcat(text, "- Starting location: ", names.v{self.loc}, "\n");

            % Movements
            text = strcat(text, "- Movements:");
            if sum(z(1:self.dv.n_e)) == 0
                text = strcat(text, " none\n");
            else
                text = strcat(text, "\n");
                for i = 1:self.dv.n_e
                    if z(i) == 1
                        [origin, destination] = findedge(Gv, i);
                        text = strcat(text, "  > From ", names.v{origin}, " to ", names.v{destination}, "\n");
                    end
                end
            end
            text = strcat(text, "- Stops at: ");
            for i = 1:self.dv.n
                if z(self.dv.n_e + i) == 1
                    text = strcat(text, names.v{i}, " ");
                end
            end
            text = strcat(text, "\n");

            % Tasks
            % TODO: This should print an ordered list of completed tasks
            text = strcat(text, "- Tasks:");
            if sum(z(self.dv.n_e + self.dv.n + (1:self.dv.n_t))) == 0
                text = strcat(text, " none\n");
            else
                text = strcat(text, "\n");
                for i = 1:self.dv.n_t
                    if z(self.dv.n_e + self.dv.n + i) == 1
                        [location, task] = find(self.idx_Evt == i);
                        text = strcat(text, "  > Task ", names.t{task}, " at location ", names.v{location}, "\n");
                    end
                end
            end

            % Location order
            % TODO: This should print an ordered list of visited locations
            text = strcat(text, "- Order (EXPERIMENTAL):");
            times = z(self.dv.n_e + self.dv.n + self.dv.n_t + (1:self.dv.n));
            [~, order] = sort(times);
            for i = 1:self.dv.n
                text = strcat(text, " ", names.v{order(i)}, ",");
            end
            text = strcat(text, "\n");

            % Add split
            text = strcat(text, "----------------\n");

        end

        function idx = get_end_node(self, z)
            % Returns the index of the graph Gv corresponding to the end-node of the agent for
            % the partial solution @z
            idx = 0;
            for i = 1:self.dv.n
                if z(self.dv.n_e + i) == 1
                    idx = i;
                    break;
                end
            end
        end

        function idx = get_movement_edges(self, z)
            % Returns the edge-indeces of the graph Gv corresponding to the movements
            % performed by the agent as indicated in the partial solution @z
            idx = find(z(1:self.dv.n_e) == 1);
        end

        function [l_idx, t_idx] = get_task_idx(self, z)
            % Returns the indeces of the location-tasks of the agent corresponding to
            % the partial solution @z
            num_tasks = sum(z(self.dv.n_e + self.dv.n + (1:self.dv.n_t)));
            if num_tasks == 0
                l_idx = [];
                t_idx = [];
            else
                l_idx = zeros(1, num_tasks);
                t_idx = zeros(1, num_tasks);
                j = 1;
                for i = 1:self.dv.n_t
                    if z(self.dv.n_e + self.dv.n + i) == 1
                        [l_idx(j), t_idx(j)] = find(self.idx_Evt == i);
                        j = j+1;
                    end
                end
            end
        end

        function info = result_info(self, z)
            % This function returns a structure containing information related to the solution @z of the agent
            % Things like an ordered list of visited locations (with timestamps), etc.

            % Tasks
            num_tasks = sum(z(self.dv.n_e + self.dv.n + (1:self.dv.n_t)));
            info.tasks.num_tasks = num_tasks;
            info.tasks.t_idx = zeros(num_tasks, 1); % Index of task in self.names.t
            info.tasks.v_idx = zeros(num_tasks, 1); % Index of tasks in self.names.v
            info.tasks.times = zeros(num_tasks, 1);

            aux_idx = 0;
            for i = 1:self.dv.n_t
                if z(self.dv.n_e + self.dv.n + i) == 1
                    aux_idx = aux_idx + 1;
                    [v_idx, t_idx] = find(self.idx_Evt == i);
                    info.tasks.v_idx(aux_idx) = v_idx;
                    info.tasks.t_idx(aux_idx) = t_idx;
                end
            end

            % Pathing and timing
            [times, locs] = sort(z(self.dv.n_e + self.dv.n + self.dv.n_t + (1:self.dv.n)));
            num_zeros = find(times >= min(self.Gv.Edges.Weight)-0.0000001, 1); % This could lead to a bug due to MILP solver tolerance
            info.times = [0; times(num_zeros:end)];
            info.v_idx = [self.loc; locs(num_zeros:end)];
            info.num_locs = length(info.times);
            
            info.tasks.Evt = self.Evt;
            info.tasks.Evt(info.tasks.Evt==1) = inf;
            task_add_time = 0;
            for l = 1:info.num_locs
                info.times(l) = info.times(l) + task_add_time;
                for i = 1:info.tasks.num_tasks
                    if info.tasks.v_idx(i) == info.v_idx(l)
                        task_time = self.param.cost_task(info.tasks.t_idx(i));
                        info.times(l) = info.times(l) + task_time;
                        task_add_time = task_add_time + task_time;
                        info.tasks.times(i) = info.times(l);
                        info.tasks.Evt(info.tasks.v_idx(i), info.tasks.t_idx(i)) = info.times(l);
                    end
                end
            end
            
        end

    end
    
end