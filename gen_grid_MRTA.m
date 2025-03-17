%% gen_grid_MRTA - Generate a random MRTA problem on a grid
%
% INPUTS:
%   @rows: Number of rows of the locations grid
%   @cols: Number of columns of the locations grid
%   @m: Number of agents
%   @k: Number of task-types
%   @opt: Structure containing several options
%         .Evt_sparse: Determines sparsity of matrix Evt. Scalar in (0, 1]
%         .Eat_sparse: Determines sparsity of matrix Eat. Scalar in (0, 1]
%   @seed: seed for RNG related to grid parameters
%
% OUTPUTS:
%   @Gv: Instance of digraph representing the location graph
%   @Evt: nxt matrix that indicates which task-types can be assigned to each location
%   @Eat: mxt matrix that indicates which agents can perform each task-type
%   @names: Structure containing the names of the locations (.v), tasks (.t) and agents (.a)
%   @param: Structure containing other information of the MRTA problem, including cost of each task
%

function [Gv, Evt, Eat, names, param] = gen_grid_MRTA(rows, cols, m, k, opt, seed)
    arguments
        rows (1, 1) {mustBeNonnegative, mustBeInteger}
        cols (1, 1) {mustBeNonnegative, mustBeInteger}
        m (1, 1) {mustBeNonnegative, mustBeInteger}
        k (1, 1) {mustBeNonnegative, mustBeInteger}
        opt struct
        seed(1, 1) {mustBeNonnegative} = 0;
    end
    if seed~=0; rng(seed); end

    %% Create the location graph
    n = rows * cols; % Number of locations

    % Location names
    names.v = cell(1, n);
    for r = 1:rows
        for c = 1:cols
            i = (r-1)*cols + c;
            names.v{i} = sprintf('v%i-%i', r, c);
        end
    end

    Ev = compute_Ev_grid(rows, cols);
    Ev = abs(10*sprandn(Ev)); % Get random weights respecting sparsity pattern
    Ev = Ev/(max(max(Ev))); % Normalize the times so that the largest one has time 1
    Gv = digraph(Ev,names.v,'omitselfloops');

    %% Define tasks and agents
    names.t = sprintfc('t%i', 1:k); % Tasks
    names.a = sprintfc('a%i', 1:m); % Agents
    
    % Evt
    Evt = rand(n, k);
    Evt(Evt >= (1- opt.Evt_sparse)) = 1;
    Evt(Evt < (1 - opt.Evt_sparse)) = 0;

    % Eat
    Eat = rand(m, k);
    Eat(Eat >= (1 - opt.Eat_sparse)) = 1;
    Eat(Eat < (1 - opt.Eat_sparse)) = 0;

    % Tasks must be doable in at least one location
    for i = 1:k
        if all(Evt(:, i) == 0)
            Evt(randi(n), i) = 1;
        end
    end
    
    % Tasks must be doable by at least one agent
    for i = 1:k
        if all(Eat(:, i) == 0)
            Eat(randi(m), i) = 1;
        end
    end

    % Every agent must be able of doing at least one task
    for i = 1:m
        if all(Eat(i, :) == 0)
            Eat(i, randi(k)) = 1;
        end
    end
    
    % Param
    param.cost_task = rand(1, k); % Random task costs

end