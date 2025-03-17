%% run_rand_grid() - Generates and solves a random grid problem
%
% INPUTS:
% @param: struct, containing the parameters for random generation
%   .row_lim: [int, int], min/max limits for number of rows of grid
%   .col_lim: [int, int], min/max limits for number of colummns of grid
%   .m_lim: [int, int], min/max limits for number of agents
%   .k_lim: [int, int], min/max limits for number of task-types
%   .Evt_sparse_lim: [float, float] in [0, 1], sparsity of Evt
%   .Eat_sparse_lim: [float, float] in [0, 1], sparsity of Eat
%   .max_tv_lim: [int, int], min/max limits for max number of each task
%   .min_tv_lim: [int, int], min/max limits for min number of each task
% @seed: struct, used to set seeds for the rangom generation
%   .graph: int, seed for generation of random graph
%   .state: int, seed for assignment of tasks and locations 
%
% NAME/VALUE INPUTS:
% 'save_me': bool (optional), determines if the results are saved in a file
% 'save_name': str (optional) no spaces, used for name generation of saved file
% 'save_to_results': bool (optional), saves the mat files to results/
% 'save_add_time': bool (otional), adds a timestamp to the name of the file
% 'run_heur': bool (optional), determines if the heuristic solver is run
% 'run_central': bool (optional), determines if the central solver is run
%
% OUTPUTS:
% @res: Structure containing results of the solver (timing results, num. iterations, exit flag, etc.)
% @test_info: Structure containing the information defining the test.
%             This information is saved to a .mat file in '/results' if save = 1
%

function [res, test_info] = run_rand_grid(param, seed, save_me, save_name, opt)
    arguments
        param struct
        seed struct
        save_me (1, 1) {mustBeNumericOrLogical} = true
        save_name {mustBeText} = ""
        opt.save_to_results = true
        opt.save_add_time = true
        opt.run_heur = true
        opt.run_central = true
        opt.ignore_seed_central = []
        opt.ignore_seed_heur = []
    end
    % Check arguments
    if strcmp(save_name, ""); save_name = "test"; end
    if any(opt.ignore_seed_heur == seed.state) && any(opt.ignore_seed_central == seed.state)
        fprintf(" Problem ignored! Graph seed: %d, State seed: %d\n\n", seed.graph, seed.state);
        return
    end

    %% Generate random grid
    rng(seed.graph); % Set seed for generating random graph
    
    % Size of grid, agents and task-types
    n_rows = randi(param.row_lim); % Rows of grid
    n_cols = randi(param.col_lim); % Columns of grid
    n = n_rows * n_cols;
    m = randi(param.m_lim); % Number of agents
    k = randi(param.k_lim); % Number of task-types

    % Sparsity of matrices Evt and Eat
    opt_MRTA.Evt_sparse = (param.Evt_sparse_lim(2) - param.Evt_sparse_lim(1))*rand(1) + param.Evt_sparse_lim(1);
    opt_MRTA.Eat_sparse = (param.Eat_sparse_lim(2) - param.Eat_sparse_lim(1))*rand(1) + param.Eat_sparse_lim(1);
    
    % Number of assigned tasks for each task-type
    opt_MRTA.max_tv = randi(param.max_tv_lim);
    opt_MRTA.min_tv = randi(param.min_tv_lim);

    % Generate MRTA problem
    [Gv, Evt, Eat, names, param_MRTA] = gen_grid_MRTA(n_rows, n_cols, m, k, opt_MRTA, seed.graph);
    [t_loc, a_loc] = gen_rand_state(Evt, names, opt_MRTA, seed.state);

    prob = MRTA(Gv, Evt, Eat, names, 'param', param_MRTA);
    num_tasks_tt = prob.get_num_tasks(t_loc);
    fprintf(" Problem with n = %d, m = %d, k = %d, num_tasks = %d, nz = %d.", prob.n, prob.m, prob.k, num_tasks_tt, prob.nz);
    fprintf(" Graph seed: %d, State seed: %d\n", seed.graph, seed.state);

    %% Save initial configuration
    save_name = sprintf("%s_gs_%s_ss_%s", save_name, num2str(seed.graph), num2str(seed.state));
    if opt.save_to_results; save_name = sprintf("results/%s", save_name); end
    if opt.save_add_time; save_name   = sprintf("%s_%s", save_name, string(datetime('now', 'Format', 'y-M-d_HH-mm-ss'))); end 

    % Generate test_info structure, which contains info related to the test setup
    test_info.type = 'grid';
    test_info.param = param;
    test_info.seed = seed;
    test_info.num_tasks = num_tasks_tt;
    test_info.n = prob.n;
    test_info.n_e = prob.n_e;
    test_info.m = prob.m;
    test_info.k = prob.k;
    test_info.nz = prob.nz;
    test_info.a_loc = a_loc;
    test_info.t_loc = t_loc;
    test_info.solved_heur = false;
    test_info.solved_central = false;

    % Save info structure
    if save_me
        save(save_name, 'test_info');
    end

    %% Solve MRTA problem using dual Lagrangian heuristic solver

    if opt.run_heur && ~any(opt.ignore_seed_heur == seed.state)

        fprintf("  > Starting dual Lagrangian heuristic solver at %s.", string(datetime('now', 'Format', 'd-MMM-y HH:mm:ss')));
        [z_d, V_d, info_d] = heur_distLag_MRTA(prob, a_loc, t_loc);
        fprintf(" Finished at %s with exit flag %d\n", string(datetime('now', 'Format', 'd-MMM-y HH:mm:ss')), info_d.e_flag);
        test_info.solved_heur = true;
    
        % Generate res structure, containing the results obtained from the solver
        res.time_d_par = info_d.run_time_par;
        res.time_d_par_total = info_d.run_time_par_total;
        res.time_d_all = info_d.run_time_all;
        res.time_d_centr = info_d.run_time_centr;
        res.e_flag_d = info_d.e_flag;
        res.feasible_d= info_d.feasible;
        res.k_d = info_d.k;
        res.V_d = V_d;
    
        % Update saved info
        if save_me
            save(save_name, '-append', 'test_info', 'res');
        end

    end

    %% Solve MRTA problem using centralized MILP solver

    if opt.run_central && ~any(opt.ignore_seed_central == seed.state)

        fprintf("  > Starting central solver at %s.", string(datetime('now', 'Format', 'd-MMM-y HH:mm:ss')));
        prob.build_MILP(a_loc, t_loc, 'sparse', true, 'condense', true); % Build MIQP problem
        [z_c, V_c, info_c] = prob.miqp.solve_mosek();
        fprintf(" Finished at %s with exit flag %d\n", string(datetime('now', 'Format', 'd-MMM-y HH:mm:ss')), info_c.e_flag);
        test_info.solved_central = true;
    
        % Update res structure with results of centralized solver
        res.time_c = info_c.run_time;
        res.e_flag_c = info_c.e_flag;
        res.feasible_c = info_c.feasible;
        res.V_c = V_c;
        
        % Update saved info
        if save_me
            save(save_name, '-append', 'test_info', 'res');
        end

    end

    %% Results if both solvers were run

    if opt.run_central && opt.run_heur

        try
            res.err_V = (V_d - V_c)/V_c;
            fprintf("  > Time centralized: %g. Time decentralized (par_total): %g. Cost centralized: %g. Cost decentralized: %g. Error: %3.2f%%\n", res.time_c, res.time_d_par_total, res.V_c, res.V_d, 100*res.err_V);
            if save_me
                save(save_name, '-append', 'test_info', 'res');
            end
        catch
            res.err_V = NaN;
            res.time_c = NaN;
            res.e_flag_c = NaN;
            res.feasible_c = NaN;
            res.V_c = NaN;
        end
    
        fprintf("\n");
    
    end

end