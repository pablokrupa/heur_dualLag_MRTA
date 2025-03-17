%% A class for defining, working and solving Multi-Robot Task-Allocation problems
%
% The class stores the problem information:
% - Location information: names, directed graph, colors, etc.
% - Agent information: number, names, position, etc.
% - Task information: types, number, locations, etc.
% - Color information for all the above (e.g., which agents can do which tasks)
%
% The class then provides methods for solving the problem and interpreting
% the results, as well as methods for updating the problem (e.g., new tasks).
% 

classdef MRTA < handle

    properties
        % For defining the problem setup
        n % Number of locations
        n_e % Number of edges of the locations graph
        nz % Number of decision variables of the MIQP problem
        m % Number of agents
        k % Number of types (colors) of tasks
        Gv % Graph that defines the locations and edges between them
        Evt % Matrix that defines the colors of each location
        Evt_full % Matrix for storing the full Evt (before any possible condensing)
        Eat % Matrix that defines the colors of each agent
        Eatv % Matrix that defines which agents can peform tasks in which locations
        names % Structure containing names of locations (v), agents (a) and tasks (t)
        task_info % Structure containing useful information about pending tasks
        
        % For defining and working with the MIQP problem
        miqp % Instance of MIQPholder that defines the MILP of the MRTA problem
        state % Structure containing the current state (agent and task locations a_loc and t_loc, etc.)
        param % Structure containing the parameters of the MIQP problem

        % Others
        opt % Options of the MRTA class itself
        agent % For storing the agents information
        idx % Contains index information of the decision variables
    end

    methods

        % Constructor
        function obj = MRTA(Gv, Evt, Eat, names, opt)
            % MRTA(Gv, Evt, Eat, names)
            % The constructor takes the graphs @Gv, @Evt and @Eat
            % that define the problem.
            % It also the structure @names, defining its names.
            arguments
                Gv digraph
                Evt (:, :)
                Eat (:, :)
                names struct
                opt.opt struct = struct();
                opt.gpad struct = struct();
                opt.param struct = struct();
            end

            obj.n = length(names.v);
            obj.n_e = size(Gv.Edges, 1);
            obj.m = length(names.a);
            obj.k = length(names.t);
            obj.Gv = Gv;
            obj.Evt = Evt;
            obj.Evt_full = Evt;
            obj.Eat = Eat;
            obj.Eatv = Eat*Evt';
            obj.nz = (obj.n_e + 2*obj.n)*obj.m + sum(sum(obj.Eatv));
            obj.names = names;
            obj.idx = struct();
            obj.state.a_loc = {};
            obj.state.t_loc = {};

            obj.opt = MRTA.default_opt();
            fn = fieldnames(obj.opt);
            % Override fields
            for i = 1:numel(fn)
                if isfield(opt.opt, fn{i})
                    obj.opt.(fn{i}) = opt.opt.(fn{i});
                end
            end

            obj.param = MRTA.default_param();
            fn = fieldnames(obj.param);
            % Override fields
            for i = 1:numel(fn)
                if isfield(opt.param, fn{i})
                    obj.param.(fn{i}) = opt.param.(fn{i});
                end
            end

        end

        function build_MILP(self, a_loc, t_loc, opt)
            % MRTA.build_MILP(a_loc, t_loc, [opt]) - Builds central MILP problem
            % 
            % Builds the central MILP problem, whose optimal solution
            % provides the optimal solution of the task-allocation problem.
            % 
            % INPUTS:
            %  @a_loc: list of length m indicating where each agent is located.
            %          Order in the list is assumed to be the same as self.names.a
            %          example: a_loc = {'v1', 'v1', 'v3'}.
            %  @t_loc: cell array of length k indicating location of tasks.
            %          Each element is itself a cell array containing location names.
            %          Order in the list is assumed to be the same as self.names.t
            %          A location can have more than one task.
            %          example (4 tasks): t_loc = {{'v1'}, {'v2', 'v4'}, {}, {'v4'}}
            %  @opt: Structure with options
            %        .sparse: Boolean that indicates if the MILP ingredients are made sparse
            %        .condense: Boolean that indicates if only decision variables for pending tasks should be added
            % 
            arguments
                self MRTA
                a_loc
                t_loc
                opt.sparse (1, 1) {mustBeNumericOrLogical} = false
                opt.condense (1, 1) {mustBeNumericOrLogical} = true
            end

            % Check arguments
            if length(a_loc) ~= self.m; error('Location of all agents must be specified'); end
            self.check_task_locations(t_loc);

            % Consense problem
            if opt.condense
                self.condense_Evt(t_loc);
            end

            % Vector of decision variables: z = [x_{i, j}^a1, t_{j, k}^a1, p_j^a1, ..., x_{i, j}^am, t_{j, k}^am, p_j^am] 

            % Construct agents
            self.agent = cell(1, self.m);
            curr_agent_idx = 0;

            for a = 1:self.m
                self.agent{a} = Agent(self.Gv, self.Evt, self.Eat(a, :), 'name', self.names.a{a}, ...
                    'loc', self.v_pos(a_loc{a}), 't_loc', t_loc, 'param', self.param, 'idx', curr_agent_idx);
                curr_agent_idx = curr_agent_idx + self.agent{a}.dv.nz; % Update position of next agent in global z       
            end

            % Construct cost function
            c_cell = cell(1, self.m);

            for a = 1:self.m
                c_cell{a} = self.agent{a}.build_c();
            end
            c_vec = cat(1, c_cell{:});

            % Equality constraints for in=out
            inout_cell = cell(1, self.m); % There is a group of constraints for each agent
            inout_vec_cell = cell(1, self.m);

            for a = 1:self.m
                [inout_cell{a}, inout_vec_cell{a}] = self.agent{a}.gen_in_out_eq_constraints();
            end

            inout_mat = blkdiag(inout_cell{:}); % Matrix that contains the equality constraint for in=out restriction
            inout_vec = cat(1, inout_vec_cell{:}); % Vector that contains the equality constraint for in=out restriction

            % Equality constraints for single-time task completion
            [task_comp_mat, task_comp_vec] = self.gen_single_task_eq_constraints(t_loc);

            % Inequality constraints for agent only performing a task if it has entered the location
            can_do_cell = cell(1, self.m); % There is a group of constraints for each agent
            can_do_vec_cell = cell(1, self.m);

            for a = 1:self.m
                [can_do_cell{a}, ~, can_do_vec_cell{a}] = self.agent{a}.gen_can_do_ineq_constraints();
            end

            can_do_mat = blkdiag(can_do_cell{:}); % Matrix that contains the inequality constraint for can_do restriction
            can_do_vec = cat(1, can_do_vec_cell{:}); % Vector that contains the inequality constraint for can_do restriction
            % Constraints are can_do_mat <= can_do_vec

            % Inequality constraints to avoid subtours by assigning a monotonicaly increasing time to each location the robot visits
            % This solution means that a robot can only visit a location (node) once
            time_cell = cell(1, self.m); % There is a group of constraints for each agent
            time_lb_cell = cell(1, self.m);
            time_ub_cell = cell(1, self.m);

            for a = 1:self.m
                [time_cell{a}, time_lb_cell{a}, time_ub_cell{a}] = self.agent{a}.gen_subtour_ineq_constraints('with_box', true); 
            end

            time_mat = blkdiag(time_cell{:}); % Matrix that contains the inequality constraint for time restriction
            time_lb_vec = cat(1, time_lb_cell{:}); % Vector containing the lower bound of the time-restriction inequalities
            time_ub_vec = cat(1, time_ub_cell{:}); % Vector containing the upper bound of the time-restriction inequalities

            % Generate matrices for forcing integer variables (for the MIQP problem)
            D_cell = cell(1, self.m);
            li_cell = cell(1, self.m);
            ui_cell = cell(1, self.m);

            for a = 1:self.m
                [D_cell{a}, li_cell{a}, ui_cell{a}] = self.agent{a}.gen_integer_constraints();
            end

            % Construct ingredients of the MIQP problem (arguments of the MIQPholder constructor)
            Q = zeros(self.nz, self.nz); % NOTE: This is a MILP problem, so Q=0
            c = c_vec;
            A = [can_do_mat; time_mat];
            ub = [can_do_vec; time_ub_vec];
            lb = [-self.n*ones(length(can_do_vec), 1);
                  time_lb_vec];
            C = [inout_mat; task_comp_mat];
            b = [inout_vec; task_comp_vec];
            D = blkdiag(D_cell{:});
            li = cat(1, li_cell{:});
            ui = cat(1, ui_cell{:});

            % Limits of integral variables
            lbz_cell = cell(1, self.m);
            ubz_cell = cell(1, self.m);
            idx_int_cell = cell(1, self.m);
            for a = 1:self.m
                [lbz_cell{a}, ubz_cell{a}] = self.agent{a}.gen_box_constraints();
                idx_int_cell{a} = self.agent{a}.get_idx_integrality();
                idx_int_cell{a} = idx_int_cell{a} + self.agent{a}.dv.idx;
            end
            lbz = cat(1, lbz_cell{:});
            ubz = cat(1, ubz_cell{:});
            idx_int = cat(2, idx_int_cell{:});

            % Build MIQP problem
            self.miqp = MIQPholder(Q, c, A, lb, ub, C, b, D, li, ui, lbz, ubz, idx_int);
            if opt.sparse
                self.miqp.to_sparse();
            end

            % Update other properties
            self.gen_task_info(t_loc); % Build task info (self.task_info)
            self.state.a_loc = a_loc;
            self.state.t_loc = t_loc;
        end

        function build_distMILP(self, a_loc, t_loc, opt)
            % MRTA.build_distMILP(a_loc, t_loc, [opt]) - Builds the distMILP problem
            % 
            % Builds an MILP problem for each agent that does not include the constraint
            % specifying that each task should only be completed once. 
            % This constraints is dealt with by the dual Lagrangian distributed algorithm.
            % 
            % INPUTS:
            %  @a_loc: list of length m indicating where each agent is located.
            %          Order in the list is assumed to be the same as self.names.a
            %          example: a_loc = {'v1', 'v1', 'v3'}.
            %  @t_loc: cell array of length k indicating location of tasks.
            %          Each element is itself a cell array containing location names.
            %          Order in the list is assumed to be the same as self.names.t
            %          A location can have more than one task.
            %          example (4 tasks): t_loc = {{'v1'}, {'v2', 'v4'}, {}, {'v4'}}
            %  @opt: Structure with options
            %        .sparse: Boolena that indicates if the MIQP ingredients are made sparse
            %        .condense: Boolean that indicates if only decision variables for pending tasks should be added
            %
            arguments
                self MRTA
                a_loc
                t_loc
                opt.param struct  = struct()
                opt.condense (1, 1) {mustBeNumericOrLogical} = true
            end
            
             % Check arguments
            if length(a_loc) ~= self.m; error('Location of all agents must be specified'); end
            self.check_task_locations(t_loc);
            if isempty(fieldnames(opt.param)); opt.param = self.param; end

            % Consense problem
            if opt.condense
                self.condense_Evt(t_loc);
            end
            
            % Construct agents
            self.agent = cell(1, self.m);
            curr_agent_idx = 0;

            for a = 1:self.m
                self.agent{a} = Agent(self.Gv, self.Evt, self.Eat(a, :), 'name', self.names.a{a}, ...
                    'loc', self.v_pos(a_loc{a}), 't_loc', t_loc, 'param', opt.param, 'idx', curr_agent_idx);
                curr_agent_idx = curr_agent_idx + self.agent{a}.dv.nz; % Update position of next agent in global z

                self.agent{a}.build_MILP();
            end

            % Update other properties
            self.gen_task_info(t_loc); % Build task info (self.task_info)
            self.state.a_loc = a_loc;
            self.state.t_loc = t_loc;

        end

        function condense_Evt(self, t_loc)
            % This method condenses the Evt matrix so that we only declare
            % decision varibles corresponding to the currently pending tasks.
            % This reduces the number of decision variables and constraints.

            self.Evt = zeros(self.n, self.k); % Initialize to zeros
            
            for i = 1:self.k
                t_loc_i = t_loc(i);
                for j = 1:length(t_loc_i{:})
                    t_pos = self.t_pos(self.names.t{i});
                    v_pos = self.v_pos(t_loc_i{1}{j});
                    if ~isempty(t_pos) && ~isempty(v_pos)
                        if self.Evt_full(v_pos, t_pos) % Assign a 1 only if in Evt_full
                            self.Evt(v_pos, t_pos) = 1;
                        end
                    end
                end
            end

            self.Eatv = self.Eat*self.Evt'; % Update Eatv matrix
            self.nz = (self.n_e + 2*self.n)*self.m + sum(sum(self.Eatv)); % Update nz

        end

        function [is_feas, feas_eq, feas_ineq, feas_int] = check_feas(self, z, tol_eq, tol_ineq)
            % Checks if a given solution @z is primal feasible
            % Feasibility of the constraints are checked using:
            %   @tol_eq for satisfaction of the equality constraints
            %   @tol_ineq for satisfaction of the inequality constraints
            arguments
                self MRTA
                z (:, 1)
                tol_eq (1, 1) {mustBeNonnegative} = 1e-2
                tol_ineq (1, 1) {mustBeNonnegative} = 1e-2
            end

            % Check equality constraints
            res_eq = self.miqp.primal.C*z - self.miqp.primal.b;
            res_ineq = self.miqp.primal.A*z;
            res_int = self.miqp.primal.D*z;

            feas_eq = norm(res_eq, inf) <= tol_eq;
            feas_ineq = all(res_ineq <= (self.miqp.primal.ub + tol_ineq)) && all(res_ineq >= (self.miqp.primal.lb - tol_ineq));
            feas_int = all(res_int == 1 | res_int == 0);
            
            is_feas = feas_eq && feas_ineq && feas_int;

        end

        function [feas_flag, list] = check_task_feas(self, z, t_loc)
            % Checks if the given solution @z satisfies all tasks listed in @t_loc being done by the agents
            % This information is provided by the output @feas_flag, which can take the following values
            %  = -1 : There are tasks in z that are not done by any robot 
            %  = -2 : There are some tasks that are not done and some that are done more than once
            %  = 0  : All tasks are done, but some tasks might be done more than once
            %  = 1  : All tasks are done exactly once
            %
            % The function also outputs list that provide information about which task-locaitons have been 
            % performed by which agents, including the number of times they have been performed
            % 

            % Get number of task-locations in t_loc
            num_tasks = self.get_num_tasks(t_loc);

            % Initialize
            z_round = round_res(self, z);
            
            % Initialize list structure
            list.t = cell(num_tasks, 1); % Task name
            list.l = cell(num_tasks, 1); % Location
            list.a = cell(num_tasks, 1); % Agents
            list.times = zeros(num_tasks, 1); % Number of times the task is performed
            
            % Retreive the list of tasks performed by each agent along with the locations
            l_idx = cell(1, self.m);
            t_idx = cell(1, self.m);

            for a = 1:self.m
                z_a = z_round(self.agent{a}.dv.idx + (1:self.agent{a}.dv.nz));
                [l_idx{a}, t_idx{a}] = self.agent{a}.get_task_idx(z_a);
            end

            % Go through list of task-locations provided by the user in t_loc
            idx_tl = 0;
            for i = 1:self.k
                t_loc_i = t_loc(i);
                t_pos = self.t_pos(self.names.t{i});
                for j = 1:length(t_loc_i{:})
                    idx_tl = idx_tl + 1;
                    list.t{idx_tl} = self.names.t{i};
                    list.l{idx_tl} = t_loc_i{1}{j};
                    l_pos = self.v_pos(t_loc_i{1}{j});
                    
                    for a = 1:self.m
                        for l = 1:length(l_idx{a})
                            if t_pos == t_idx{a}(l) && l_pos == l_idx{a}(l)
                                list.times(idx_tl) = list.times(idx_tl) + 1;
                                if list.times(idx_tl) == 1
                                    list.a{idx_tl} = self.names.a{a};
                                else
                                    list.a{idx_tl} = strcat(list.a{idx_tl}, {', '}, self.names.a{a});
                                end
                            end
                        end
                    end

                end
            end

            % Determine value of output flag
            if all(list.times == 1)
                feas_flag = 1;
            elseif any(list.times == 0)
                feas_flag = -1;
            else
                feas_flag = 0;
            end

            if feas_flag == -1 && max(list.times) > 1
                feas_flag = -2;
            end

        end

        function passed = check_task_locations(self, t_loc, opt)
            % MRTA.check_tast_locations(t_loc) - Checks if the tasks can be assigned
            % to the given locations
            arguments
                self MRTA
                t_loc
                opt.verbose {mustBeNumericOrLogical} = self.opt.verbose
            end
            
            passed = true;
            for i = 1:self.k
                t_loc_i = t_loc(i);
                for j = 1:length(t_loc_i{:})
                    if ~self.check_task_location(self.names.t{i}, t_loc_i{1}{j}, 'verbose', opt.verbose)
                        passed = false;
                        break;
                    end
                end
                if ~passed; break; end
            end

        end

        function passed = check_task_location(self, t, v, opt)
            % MRTA.check_tast_location(t, v, ['verbose'])
            % Checks if the location name @v has task name @t
            arguments
                self MRTA
                t string
                v string
                opt.verbose {mustBeNumericOrLogical} = self.opt.verbose
            end

            t_pos = self.t_pos(t);
            v_pos = self.v_pos(v);
            if isempty(t_pos) || isempty(v_pos)
                passed = false; % Location @v or task @t does not exist
            else
                passed = self.Evt(v_pos, t_pos);
            end
            if opt.verbose && ~passed
                warning(append('Task ', t, ' cannot be assigned to location ', v));
            end

        end

        function num_tasks = get_num_tasks(self, t_loc)
            % Returns the number of tasks in t_loc
            num_tasks = 0;
            for i = 1:self.k
                t_loc_i = t_loc(i);
                num_tasks = num_tasks + length(t_loc_i{:});
            end
        end

        function gen_task_info(self, t_loc)
            % Generates a structure containing useful information about the pending tasks
            % Saves this information in property sef.task_info
            
            num_tasks = self.get_num_tasks(t_loc);
            self.task_info.num_tasks = num_tasks; % Number of pending tasks

            % List of pending tasks and their respective locations
            self.task_info.t = cell(num_tasks, 1); % Name of tasks
            self.task_info.v = cell(num_tasks, 1); % Name of locations
            self.task_info.t_idx = cell(num_tasks, 1); % Index of task in self.names.t
            self.task_info.v_idx = cell(num_tasks, 1); % Index of tasks in self.names.v

            aux_idx = 0;
            for i = 1:self.k
                t_loc_i = t_loc(i);
                t_pos = self.t_pos(self.names.t{i});
                for j = 1:length(t_loc_i{:})
                    aux_idx = aux_idx + 1;
                    self.task_info.t{aux_idx} = self.names.t{i};
                    self.task_info.t_idx{aux_idx} = t_pos;
                    self.task_info.v{aux_idx} = t_loc_i{1}{j};
                    self.task_info.v_idx{aux_idx} = self.v_pos(t_loc_i{1}{j});
                end
            end

            % Agent information
            self.task_info.num_agents = zeros(num_tasks, 1); % Stores how many agents can do the task
            self.task_info.a = cell(num_tasks, 1); % Stores name of agents that can do the task
            self.task_info.a_idx = cell(num_tasks, 1); % Stores name of agents that can do the task

            Evt_all = zeros(self.n, self.k);
            for a = 1:self.m
                Evt_all = Evt_all + self.agent{a}.Evt;
            end
            for i = 1:num_tasks
                    self.task_info.num_agents(i) = Evt_all(self.task_info.v_idx{i}, self.task_info.t_idx{i});
            end
            for a = 1:self.m
                for i = 1:num_tasks
                    if self.agent{a}.Evt(self.task_info.v_idx{i}, self.task_info.t_idx{i}) == 1
                        self.task_info.a{i}{end+1} = self.agent{a}.name;
                        self.task_info.a_idx{i}{end+1} = a;
                    end
                end
            end

        end

        function idx = get_task_location_idx(self, t, v, opt)
            % Returns the index of the task-location @t-@v in the
            % equalities for single-time task completion
            arguments
                self MRTA
                t string
                v string
                opt.verbose {mustBeNumericOrLogical} = false
            end
            if self.check_task_location(t, v, 'verbose', opt.verbose)
                t_pos = self.t_pos(t);
                v_pos = self.v_pos(v);
                idx = self.idx.Evt(v_pos, t_pos);
            else
                idx = [];
            end
            
        end

        function passed = check_agent_locations(self, a_loc, opt)
            arguments
                self MRTA
                a_loc
                opt.verbose {mustBeNumericOrLogical} = self.opt.verbose
            end
            
            passed = true;
            for i = 1:self.m
                if ~self.check_agent_location(self.names.a{i}, a_loc{i}, 'verbose', opt.verbose)
                    passed = false;
                    break;
                end
            end

        end

        function passed = check_agent_location(self, a, v, opt)
            % MRTA.check_agent_location(a, v, ['verbose'])
            % Checks if the agent name @a can be assigned to location name @v
            arguments
                self MRTA
                a string
                v string
                opt.verbose {mustBeNumericOrLogical} = self.opt.verbose
            end

            a_pos = self.a_pos(a);
            v_pos = self.v_pos(v);
            if isempty(a_pos) || isempty(v_pos)
                passed = false; % Location @v or agent @a does not exist
            else
                passed = true;
            end
            if opt.verbose && ~passed
                warning(append('Agent ', a, ' cannot be assigned to location ', v));
            end

        end

        function pos = v_pos(self, v)
            % Find possition of name @v in list of location names
            pos = find(strcmp(self.names.v, v));
        end

        function pos = a_pos(self, a)
            % Find possition of name @a in list of agent names
            pos = find(strcmp(self.names.a, a));
        end

        function pos = t_pos(self, t)
            % Find possition of name @t in list of task names
            pos = find(strcmp(self.names.t, t));
        end

        function [loc, loc_num] = a_loc(self, a)
            % Returns the current location of the agent @a
            if ~isnumeric(a); a = self.a_pos(a); end
            loc_num = self.agent{a}.loc;
            loc = self.names.v{loc_num};
        end

        function z_a = get_z_agent(self, z, a)
            % Returns the portion of @z coresponding to agent number @a
            z_a = z(self.agent{a}.dv.idx + (1:self.agent{a}.dv.nz));          
        end

        function z_round = round_res(self, z)
            % This function returns a integral-rounded variable z_round from an approximate z
            % That is, the numbers of z corresponding to integral variables are set to 0 or 1
            z_round = zeros(length(z), 1);
            for a = 1:self.m
                z_a = get_z_agent(self, z, a);
                z_round(self.agent{a}.dv.idx + (1:self.agent{a}.dv.nz)) = self.agent{a}.round_z(z_a); 
            end
        end

        function res2tex(self, z)
            % This function prints text that interprets the decision variables @z
            z_round = round_res(self, z); % Round solution so that we have integers
            fprintf("----------------\n         Results\n----------------\n")
            for a = 1:self.m
                interpret_MIQP_agent_result(self, z_round, a);
            end
        end

        function state2tex(self)
            % This function prints text that interprets the current state (agent and task locations)
            text = "----------------\n         State\n----------------\n";

            % Agent locations
            text = strcat(text, "Agent locations: \n");
            for a = 1:self.m
                text = strcat(text, "  > Agent ", self.agent{a}.name, " at ", self.state.a_loc{a}, "\n");
            end

            % Task locations
            text = strcat(text, "Task locations: \n");
            for i = 1:self.k
                text = strcat(text, "  > Task ", self.names.t{i}, ":");
                tasks_i = self.state.t_loc{i};
                if isempty(tasks_i)
                    text = strcat(text, " none\n");
                else
                    for j = 1:length(tasks_i)
                        text = strcat(text, " ", tasks_i{j}, ",");
                    end
                    text = strcat(text, "\n");
                end
            end

            fprintf(text);
        end

        function info2tex(self, info)
            % This function prints the information about the solver information in hunam-readable form
            text = "----------------\n      Solver info\n----------------\n";
            text = strcat(text, "- Exit status: ");
            switch info.e_flag
                case 1
                    text = strcat(text, "Optimal solution found (e_flag = 1)\n");
                case -1
                    text = strcat(text, "Infeasible problem detected (e_flag = -1)\n");
                case -2
                    text = strcat(text, "Maximum number of iterations (e_flag = -2)\n");
                otherwise
                    text = strcat(text, "ERROR: Unknown value of e_flag\n");
            end
            
            text = strcat(text, "- Run-time (sec): ", num2str(info.run_time), "\n");
            text = strcat(text, "- Decision varaibles: ", num2str(self.nz), "\n");
            text = strcat(text, "- Outer iterations: ", num2str(info.iters), "\n");
            text = strcat(text, "- Inner iterations: ", num2str(info.inner_iters), "\n");
            text = strcat(text, "- Optimal value: ", num2str(info.V_opt), "\n");

            fprintf(text);

        end

        function interpret_MIQP_agent_result(self, z, a_num)
            % This function prints the text corresponding to agent @a_num for the decision variables @z
            z_a = z(self.agent{a_num}.dv.idx + (1:self.agent{a_num}.dv.nz));
            text = self.agent{a_num}.res2tex(z_a, self.Gv, self.names);
            fprintf(text);   
        end

        function ph = plot_Gv(self, opt)
            % This function plots the locations graph Gv in figure @fig_num
            % It returns a handler to the figure
            arguments
                self MRTA
                opt.fig_num (1, 1) = 1 % Figure number
                opt.plot_weights (1, 1) {mustBeNumericOrLogical} = true % Plots the weight in the graph edges
            end
            
            figure(opt.fig_num); clf(opt.fig_num);
            if opt.plot_weights
                ph = plot(self.Gv, 'EdgeLabel', self.Gv.Edges.Weight, 'Layout', 'force'); 
            else
                ph = plot(self.Gv, 'Layout', 'force'); 
            end
            ph.LineWidth = 1.0;
            ph.LineStyle = '-.';
            ph.Marker = 's';
            ph.MarkerSize = 6;
            ph.ArrowPosition = 0.9;
            ph.ArrowSize = 9;
            ph.EdgeFontSize = 11;
            ph.NodeFontSize = 13;
        end

        function plot_result(self, z, opt)
            % Plots the results for each agent from the solution vector @z
            % Each agent is plotted in a separate figure
            % The figure plots:
            %   - Starting location: green node with (s) added to the node label
            %   - End location: red node with (e) added to the node label
            %   - Tasks completed as location: with [task_name] added to the node label
            %   - Travelled edges: Marked in dark red
            arguments
                self MRTA
                z (:, 1)
                opt.plot_weights (1, 1) {mustBeNumericOrLogical} = true % Plots the weight in the graph edges
                opt.fig_num = 0
                opt.max_fig = inf;
            end

            z_round = round_res(self, z); % Round solution so that we have integers

            for a = 1 : min(self.m, opt.max_fig)

                z_a = z_round(self.agent{a}.dv.idx + (1:self.agent{a}.dv.nz));

                % Get initial basic plot of the layout
                ph = self.plot_Gv('fig_num', a + opt.fig_num, 'plot_weights', opt.plot_weights);
                
                % Highlight completed tasts
                [l_idx, t_idx] = self.agent{a}.get_task_idx(z_a);
                num_tasks = length(l_idx);
                for i = 1:num_tasks
                    try
                        loc_label = ph.NodeLabel{l_idx(i)};
                    catch ME
                        loc_label = '';
                    end
                    labelnode(ph, l_idx(i), char(strcat(loc_label, " [", self.names.t(t_idx(i)), "]")));
                end

                % Highlight starting point
                highlight(ph, self.state.a_loc{a}, 'NodeColor', 'g');
                start_label = ph.NodeLabel{self.v_pos(self.state.a_loc{a})};
                ph.NodeLabel{self.v_pos(self.state.a_loc{a})} = [start_label, ' (s)'];

                % Highlight end point
                end_node_idx = self.agent{a}.get_end_node(z_a);
                if end_node_idx > 0
                    end_label = ph.NodeLabel{end_node_idx};
                    ph.NodeLabel{end_node_idx} = [end_label, ' (e)'];
                    highlight(ph, end_node_idx, 'NodeColor', 'r');
                end

                % Highlight edges where agent moves through
                move_edges = self.agent{a}.get_movement_edges(z_a);
                if ~isempty(move_edges)
                    highlight(ph, 'Edges', move_edges', 'EdgeColor', [0.6350 0.0780 0.1840], 'LineStyle', '-', 'LineWidth', 1.4);
                end

            end

        end

        function [Ct, dt, idx] = gen_single_task_eq_constraints(self, t_loc)
            % This function returns matrix Ct and vector dt for the equality
            % constraints that force a task to be completed by a single agent
            n_t = sum(sum(self.Evt));
            task_comp_cell = cell(1, n_t); % There is a group of constraints for each location-task
            idx.t = zeros(1, n_t);
            idx.v = zeros(1, n_t);
            
            tv = 1; % Counter for the current equality index
            for t = 1:self.k
                for v = 1:self.n
                    if self.Evt(v, t) ~= 0
                        task_comp_cell{tv} = zeros(1, self.nz);
                        for a = 1:self.m
                            [has_task, idx_tv] = self.agent{a}.has_task_loc(t, v);
                            if has_task
                                task_comp_cell{tv}(idx_tv) = 1;
                            end
                        end
                        % Here we save the information of which row of the equality-constraint matrix corresponds to this
                        % task-location pair. We store this in self.idx.Evt, which mimics the Evt matrix
                        self.idx.Evt(v, t) = tv;
                        idx.t(tv) = t;
                        idx.v(tv) = v;
                        tv = tv + 1;
                    end
                end
            end

            % Matrix that contains the equality constraints for single-time task completion
            Ct = reshape([task_comp_cell{:}], self.nz, sum(sum(self.Evt)))';

            % Vector that contains the equality constraints for single-time task completion
            dt = zeros(sum(sum(self.Evt)), 1);
            for i = 1:self.k
                t_loc_i = t_loc(i);
                for j = 1:length(t_loc_i{:})
                    idx_t = self.get_task_location_idx(self.names.t{i}, t_loc_i{1}{j});
                    if ~isempty(idx_t)
                        dt(idx_t) = 1;
                    end
                end
            end

        end

    end % End of public method

    methods (Static)
        
        function opt = default_opt()
            % Returns the default options structure of the MRTA class
            opt.verbose = false;
        end

        function param = default_param()
            % Returns the default parameters of the MRTA optimization problem
            param.cost_stop = 1e-2; % Cost assigned to stopping at any location
            param.cost_task = []; % Array that indicated the cost (time) of each task
            param.cost_times = 1e-5; % Cost assigned to the decision variables that establish time
        end

    end % End of Static methods

end