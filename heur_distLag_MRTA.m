%% A heuristic distributed dual Lagrangian method for MRTA problems
%
% For intormation on the base method please see:
% Vujanic, Robin, et al. "Large scale mixed-integer optimization: A solution
% method with supply chain applications." 22nd Mediterranean Conference on 
% Control and Automation. IEEE, 2014.
%
% The base method is modifies with heuristics that are designed to
% improve the practical convergence of the algorithm for MRTA problems.
%
% The dual Lagrandian method divides the central MILP problem into small
% scale MILP problems, one for each agent. Dual variables then move the
% indvidual problems towards a consensus.
% The heuristic scheme modifies the dual variables of each agent to improve
% the overall performance of the algorithm.
% 
% Each agent's individual MILP problem is solved using a comercial MILP solver.
%
% INPUTS:
%   @mrta: Instance of MRTA
%   @a_loc: Cell array that indicates the current location of each agent
%   @t_loc: Cell array of cell arrays that indicates which locations have an assigned task
%   @opt: Structure containing the options of the solver:
%       .cost_task_mag: Number that multiplies param.cost_tasks in each agent's MILP problem
%       .alpha: Value used to compute the step-size of the dual Lagrangian method
%       .beta: Hyperparameter of the dual update heuristics
%       .fixed_s_iters: Number of iterations where s_k = alpha
%       .k_max: Maximum number of iterations
%       .k_early: Number of iterations after which we allow early termination
%       .k_heur: Number of iterations after which we start applying heuristics
%       .tol: Exit tolerance for the distributed algorithm
%       .mosek: Options structure for MOSEK solver
%   @param: Structure containing parameters of the MRTA problem
%
% OUTPUTS:
%   @z_out: Solution returned from the heuristic dual Lagrangian method
%   @V_out: Value of the objective function for @z_out
%   @into: Structure containing additional information (computation times, dual variables, number of iterations, etc.)
%

function [z_out, V_out, info] = heur_distLag_MRTA(mrta, a_loc, t_loc, opt, param)
    arguments
        mrta MRTA
        a_loc
        t_loc
        opt struct = struct()
        param struct = struct()
    end

    % Parse options
    par = inputParser;
    par.addParameter('cost_task_mag', -0.5*sum(mrta.Gv.Edges.Weight)/sqrt(mrta.n)); % Number that multiplies param.cost_tasks in each agent's MILP problem
    par.addParameter('alpha', 2*sum(mrta.param.cost_task)*mrta.m); % Value used to compute the step-size of the dual Lagrangian method
    par.addParameter('beta', 1/mrta.m); % Hyperparameter of the dual update heuristics
    par.addParameter('fixed_s_iters', 2); % Number of iterations where s_k = alpha
    par.addParameter('k_max', 50); % Maximum number of iterations
    par.addParameter('k_early', 20); % Number of iterations after which we allow early termination
    par.addParameter('k_heur', 5); % Number of iterations after which we start applying heuristics
    par.addParameter('tol', 0.001); % Exit tolerance for the distributed algorithm
    par.addParameter('mosek', struct()); % Options structure for MOSEK solver
    par.addParameter('condense', true); % Determines if the optimization problem is condensed
    par.addParameter('warmstart', false); % Apply warmstart to mosek
    
    par.StructExpand=true;
    parse(par, opt);
    opt = par.Results;

    %% Construct MILP problem of each agent
    if isempty(fieldnames(param)); param = mrta.param; end

    mrta.build_distMILP(a_loc, t_loc, 'param', param, 'condense', opt.condense);

    % Build the inequality constraint that couples the MILP problems of the agents
    [H, b, idx_H] = mrta.gen_single_task_eq_constraints(t_loc);
    % Normalize H
    H = H/mrta.m;
    b = b/mrta.m;

    %% Initialize subgradient method
    k = 0; % Iteration counter
    n_d = size(H, 1); % Number of dual variables
    lambda = zeros(n_d, 1); % Initialize dual variables lambda
    res_k = zeros(n_d, 1);
    res_k1 = res_k;
    done = false; % Used to determine when to exit the iterative algorithm
    found_semifeas = false; % Used to determine if we have found a semi-feasible solution (all tasks as done)
    info.run_time_par = 0;
    info.run_time_sum = 0;
    info.run_time_a = zeros(1, mrta.m);
    info.max_times = zeros(1, mrta.m);
    info.avrg_times_a = zeros(1, mrta.m);
    info.hist_run_times_a = cell(1, mrta.m);
    z_a = cell(1, mrta.m);
    V_a = cell(1, mrta.m);
    info_a = cell(1, mrta.m);
    for a = 1:mrta.m
        z_a{a} = NaN(mrta.agent{a}.dv.nz, 1);
        info.hist_run_times_a{a} = zeros(1, opt.k_max);  
    end
    z_best = NaN(mrta.nz, 1);
    lambda_best = NaN(n_d, 1);
    V_best = Inf;

    % Warmstart
    warmstart = cell(1, mrta.m);
    for i = 1:mrta.m; warmstart{i} = struct(); end

    % Store original c vectors for each agent and separate the H part of each agent
    c_a = cell(1, mrta.m);
    H_a = cell(1, mrta.m);
    for a = 1:mrta.m
        c_a{a} = mrta.agent{a}.miqp.primal.c;
        H_a{a} = sparse(H(:, mrta.agent{a}.dv.idx + (1:mrta.agent{a}.dv.nz)));
    end
    H = sparse(H);
    c = cat(1, c_a{:});
    % Update the c_a vectors with cost_task_mag
    for a = 1:mrta.m
        c_a{a}(mrta.agent{a}.dv.n_e + mrta.agent{a}.dv.n + (1:mrta.agent{a}.dv.n_t)) = opt.cost_task_mag*c_a{a}(mrta.agent{a}.dv.n_e + mrta.agent{a}.dv.n + (1:mrta.agent{a}.dv.n_t));
    end

    % Matrix containing best times for each task
    Evt_best = mrta.Evt;
    Evt_best(Evt_best == 1) = inf;
    Evt_a = cell(1, mrta.m);
    Evt_a_best = cell(1, mrta.m);
    lambda_a = cell(1, mrta.m);
    inc_sk = cell(1, mrta.m);
    for a = 1:mrta.m
        Evt_a{a} = mrta.agent{a}.Evt;
        Evt_a{a}(Evt_a{a} == 1) = inf;
        Evt_a_best{a} = Evt_a{a};
        lambda_a{a} = zeros(n_d, 1);
        inc_sk{a} = zeros(n_d, 1);
    end
    
    %% Run algorithm
    start_MILP_dist = tic;
    info.run_time_centr = 0; % For storing the time required for the centralized operations
    while ~done
        k = k+1;
        
        % Solve the MILP problem of each agent
        for a = 1:mrta.m
            
            lambda_inc = H_a{a}'*lambda_a{a};
            mrta.agent{a}.miqp.primal.c = c_a{a} + lambda_inc; % Update c vector of the agent's MILP problem
            [z_a{a}, V_a{a}, info_a{a}] = mrta.agent{a}.miqp.solve_mosek(opt.mosek, warmstart{a}); % Solve MILP problem using MOSEK
            info.run_time_a(a) = info.run_time_a(a) + info_a{a}.run_time;
            info.hist_run_times_a{a}(k) = info_a{a}.run_time;
            z_a{a} = mrta.agent{a}.round_z(z_a{a});
            time_agent_a = tic;
            
            % Warmstart for next iteration
            if opt.warmstart 
                warmstart{a} = info_a{a}.mosek.sol.int;
            end
            
            % Get timing information
            info_temp = mrta.agent{a}.result_info(z_a{a});
            Evt_a{a} = info_temp.tasks.Evt;

            % Update best times
            for i = 1:mrta.agent{a}.dv.n_t
                idx_v = mrta.agent{a}.task_idx.v(i);
                idx_t = mrta.agent{a}.task_idx.t(i);
                Evt_best(idx_v, idx_t) = min(Evt_best(idx_v, idx_t), Evt_a{a}(idx_v, idx_t));
                Evt_a_best{a}(idx_v, idx_t) = min(Evt_a_best{a}(idx_v, idx_t), Evt_a{a}(idx_v, idx_t));
            end

            info.run_time_centr = info.run_time_centr + toc(time_agent_a);
        end
        time_coordinator = tic;
        z_k = cat(1, z_a{:}); % Concatenate solution of all agents
        
        % Compute residual
        res_k = H*z_k - b;
        if norm(res_k - res_k1) < opt.tol
            opt.beta = 2*opt.beta;
        end
        res_k1 = res_k;

        % Compute dual varibales
        if k <= opt.fixed_s_iters
            s_k = opt.alpha;
        else
            s_k = opt.alpha/(k - opt.fixed_s_iters);
        end
        % lambda = lambda + s_k*res_k; % This would be without heuristics
        
        % Update lambda (using inc_sk if k >= opt.k_heur)
        if k < opt.k_heur
            for a = 1:mrta.m
                lambda_a{a} = lambda_a{a} + s_k*res_k;
            end
        else
            for a = 1:mrta.m
                inc_sk{a} = zeros(n_d, 1);
            end
            for i = 1:n_d

                % Go through each one of the pending tasks and act in relation to how many agents want to do it
                if res_k(i) > 0 % Task done by more than one agent

                    % Compute minimum and maximum time of each agent (between agents that want to do it)
                    t_min = inf;
                    t_max = 0;
                    a_min = 0;
                    for a = 1:mrta.m
                        time_a = Evt_a{a}(idx_H.v(i), idx_H.t(i));
                        if time_a > 0 && time_a < inf
                            if norm(time_a - t_min) < 0.001
                                time_a = 2*time_a; % To break symmetries
                            end
                            t_max = max(t_max, time_a);
                            if time_a <= t_min
                                t_min = time_a;
                                a_min = a;
                            end
                        end
                    end

                    % Give agents that want to do a task a positive or negative boost depending on their time vs. minimum time
                    for a = 1:mrta.m
                        if Evt_a{a}(idx_H.v(i), idx_H.t(i)) > 0 && Evt_a{a}(idx_H.v(i), idx_H.t(i)) < inf
                            
                            inc_sk{a}(i) = opt.beta*(Evt_a{a}(idx_H.v(i), idx_H.t(i)) - t_min)/t_max;
                            
                            % If the agent is the one with minimum time give him a boost in wanting to do that task
                            if a == a_min
                                inc_sk{a}(i) = -opt.beta*res_k(i)*mrta.m;
                            end
                        end
                    end

                elseif res_k(i) < 0 % No agent wants to do the task

                    % Compute minimum and maximum time of each agent (from past iterations)
                    t_min = inf;
                    t_max = 0;
                    a_min = 0;
                    for a = 1:mrta.m
                        time_a = Evt_a_best{a}(idx_H.v(i), idx_H.t(i));
                        if time_a > 0 && time_a < inf
                            t_max = max(t_max, time_a);
                            if time_a <= t_min
                                t_min = time_a;
                                a_min = a;
                            end
                        end
                    end

                    for a = 1:mrta.m
                        % Give agents that can do the task a random boost to wanting to do it (the a_min breaks the symmetry)
                        if Evt_a{a}(idx_H.v(i), idx_H.t(i)) == inf
                            if a == a_min
                                inc_sk{a}(i) = opt.beta*res_k(i)*mrta.m;
                            else
                                inc_sk{a}(i) = -opt.beta*res_k(i);
                            end
                        end
                    end

                else % Task is being done by exactly one agent
                    for a = 1:mrta.m

                        % Give agent a doing the task a boost
                        if Evt_a{a}(idx_H.v(i), idx_H.t(i)) > 0 && Evt_a{a}(idx_H.v(i), idx_H.t(i)) < inf
                            inc_sk{a}(i) = -1/opt.beta*res_k(i);

                        % Give agents not doing the task a penalty
                        elseif Evt_a{a}(idx_H.v(i), idx_H.t(i)) == inf
                            inc_sk{a}(i) = 1/opt.beta*res_k(i);
                        end

                    end
                end
            end
            for a = 1:mrta.m
                lambda_a{a} = lambda_a{a} + s_k*(res_k + inc_sk{a});
            end
        end

        % Store best semi-feasible solution
        [feas_flag, task_info] = mrta.check_task_feas(z_k, t_loc);
        if feas_flag == 0 && k > 1
            found_semifeas = true;
            V_k = c'*z_k;
            if V_k < V_best
                z_best = z_k;
                V_best = V_k;
                lambda_best = lambda;
                info.k_best = k;
            end
        end
        
        % Compute exit condition
        if norm(res_k, inf) <= opt.tol
            done = true;
            info.e_flag = 1;
            info.feasible = true;
        elseif k >= opt.k_max
            done = true;
            info.e_flag = -1;
            info.feasible = false;
        elseif k >= opt.k_early && found_semifeas
            done = true;
            info.e_flag = 2;
            info.feasible = false;
            z_k = z_best;
            lambda = lambda_best;
        end

        info.run_time_centr = info.run_time_centr + toc(time_coordinator);

    end
    info.run_time_all = toc(start_MILP_dist);

    %% Return solution
    if info.e_flag == -1 && V_best < Inf
        info.e_flag = 0;
        z_k = z_best;
        lambda = lambda_best;
    end

    z_out = mrta.round_res(z_k);
    V_out = c'*z_out;

    if isfield(mrta.miqp, 'primal')
        [info.feasible, feas_eq, feas_ineq, feas_int] = mrta.check_feas(z_out);
        info.feas_report.feas_eq = feas_eq;
        info.feas_report.feas_ineq = feas_ineq;
        info.feas_report.feas_int = feas_int;
    end
    info.lambda = lambda;
    info.k = k;
    for a = 1:mrta.m
        info.hist_run_times_a{a} = info.hist_run_times_a{a}(1:k);
        info.max_times(a) = max(info.hist_run_times_a{a});
        info.avrg_times_a(a) = info.run_time_a(a)/k;
    end
    info.run_time_par = max(info.run_time_a);
    info.run_time_par_total = info.run_time_par + info.run_time_centr;
    info.Evt_best = Evt_best;

end