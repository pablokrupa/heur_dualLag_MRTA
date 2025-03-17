%% A class for defining, solving and working with MIQP problems
% 
% The primal MIQP optimization problem is given by
% 
% min_z 0.5*z'*Q*z + c'*z
%  s.t. lb <= A*z <= ub
%  lbz <= z <= ubz
%  C*z = b
%  D_i*z \in {li_i, ui_i}, i = 1,...,p
%
% where dim(z) = n, dim(A) = [m, n], dim(C) = [q, n], dim(D) = [p, n].
%

classdef MIQPholder < handle

    properties
        primal % For storing the primal ingredients
        n % Dimension of primal decision variables
        m % Dimension of primal inequality constraints
        q % Dimension of primal equality constraints
        p % Dimension of number of integer constraints
        nd % Dimension of dual variables
    end

    methods
        
        % Constructor
        function obj = MIQPholder(Q, c, A, lb, ub, C, b, D, li, ui, lbz, ubz, idx_int)
            % MIQPholder(Q, c, A, lb, ub, C, b, D, li, ui)
            % The constructor simply takes the ingredients of the MIQP problem
            arguments
                Q (:, :)
                c (:, 1)
                A (:, :)
                lb (:, 1) = []
                ub (:, 1) = []
                C (:, :) = []
                b (:, 1) = []
                D (:, :) = []
                li (:, 1) = []
                ui (:, 1) = []
                lbz (:, 1) =[]
                ubz (:, 1) = []
                idx_int (:, 1) = []
            end

            % Get sizes
            obj.n = size(Q, 1);
            obj.m = size(A, 1);
            obj.q = size(C, 1);
            obj.p = size(D, 1);
            obj.nd = 2*obj.m + 2*obj.p + obj.q;

            % Give values to empty vectors
            if isempty(lb); lb = -inf*ones(obj.m, 1); end
            if isempty(ub); ub = inf*ones(obj.m, 1); end
            if isempty(b); b = zeros(obj.q, 1); end
            if isempty(li); li = zeros(obj.p, 1); end
            if isempty(ui); ui = ones(obj.p, 1); end
            if isempty(lbz); lbz = -inf*ones(obj.n, 1); end
            if isempty(ubz); ubz = inf*ones(obj.n, 1); end
            if isempty(idx_int)
                idx_non_int = find(all(D==0));
                idx_int = setdiff(1:obj.n, idx_non_int);
            end
            
            % Asign variables
            obj.primal.Q = Q;
            obj.primal.c = c;
            obj.primal.A = A;
            obj.primal.lb = lb;
            obj.primal.ub = ub;
            obj.primal.C = C;
            obj.primal.b = b;
            obj.primal.D = D;
            obj.primal.li = li;
            obj.primal.ui = ui;
            obj.primal.lbz = lbz;
            obj.primal.ubz = ubz;
            obj.primal.idx_int = idx_int;
            
        end

        function [z_out, V_out, info] = solve_mosek(self, opt, warmstart)
            % Method for solving the MILP problem using MOSEK
            % Note that this method assumes that the object defines a MILP
            % problem, i.e., that the Q matrix is zero.
            arguments
                self MIQPholder
                opt struct = struct() % TODO: Add options for MOSEK solver
                warmstart struct = struct()
            end

            if any(self.primal.Q)
                error("MIQPholder:not_MILP", "MIQPholder.solve_mosek() can only solve MILP problems");
            end

            % Parse options
            par = inputParser;
            par.addParameter('MSK_IPAR_LOG', 0); % Verbose level for MOSEK
            
            par.StructExpand=true;
            if ~isfield(opt, 'mosek'); opt.mosek = struct(); end
            parse(par, opt.mosek);
            opt.mosek = par.Results;

            % Construcy prob structure for MOSEK
            prob.c = self.primal.c;
            prob.a = [self.primal.C; self.primal.A];
            prob.blc = [self.primal.b; self.primal.lb];
            prob.buc = [self.primal.b; self.primal.ub];
            prob.blx = self.primal.lbz;
            prob.bux = self.primal.ubz;
            prob.ints.sub = self.primal.idx_int;

            % Warmstart
            fn = fieldnames(warmstart);
            for i = 1:numel(fn)
                prob.sol.int.(fn{i}) = warmstart.(fn{i});
            end
            
            % Call MOSEK solver and return solutionsp
            if opt.mosek.MSK_IPAR_LOG > 0
                [r, res] = mosekopt('minimize info', prob, opt.mosek);
            else
                [r, res] = mosekopt('minimize info echo(0)', prob, opt.mosek);
            end
            info.mosek_e_flag = r;
        
            if r == 0
                info.e_flag = 1; % Optimal solution found
                info.feasible = 1;
                z_out = res.sol.int.xx;
                V_out = self.primal.c'*z_out;
            elseif r > 0
                z_out = NaN(self.n, 1);
                V_out = NaN;
                info.feasible = 0;
                info.e_flag = -1; % Maximum iterations or time
            else
                error("Unknown MOSEK exit flag");
            end
            
            info.run_time = res.info.MSK_DINF_OPTIMIZER_TIME;
            info.mosek.info = res.info;
            info.mosek.sol.int = res.sol.int;
            info.mosek.rcodestr = res.rcodestr;
            if isfield(info.mosek, 'solsta')
                info.mosek.prosta = res.sol.int.prosta;
                info.mosek.solsta = res.sol.int.solsta;
            else
                info.mosek.prosta = 'ERROR';
                info.mosek.solsta = 'ERROR';
            end

        end

        function to_sparse(self)
            % Converts the main ingredients of the MIQP problem to sparse matrices
            self.primal.A = sparse(self.primal.A);
            self.primal.C = sparse(self.primal.C);
            self.primal.D = sparse(self.primal.D);
        end

        function int_feas = is_integer_feasible(self, z, tol)
            % is_integer_feasible(z, [tol])
            % Returns true if @z is integer feasible for the given @tol
            arguments
                self MIQPholder
                z (:, 1)
                tol double = 1e-6
            end
            int_feas = true;
            t = self.primal.D * z;
            for i = 1:self.p
                if abs(t(i) - self.primal.li(i)) > tol && abs(t(i) - self.primal.ui(i)) > tol
                    int_feas = false;
                    break;
                end
            end
        end

    end % End of public methods

end