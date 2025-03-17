%% gen_rand_state - Generate a random state for a MRTA problem
% 
% The state is the location of each agent and the pending tasks.
%
% INPUTS:
%   @Evt: nxt matrix that indicates which task-types can be assigned to each location
%   @names: Structure containing the names of the locations (.v), tasks (.t) and agents (.a)
%   @opt: Structure containing several options
%         .max_tv: Maximum number of times each task-type can be assigned overall. Positive integer
%         .min_tv: Minimum number of times each task-type can be assigned overall. Positive integer <= max_tv
%   @seed: seed for RNG related to agent and task locations
%
% OUTPUTS:
%   @t_loc: Cell array of cell arrays that indicates which locations have an assigned task
%   @a_loc: Cell array that indicates the current location of each agent
%

function [t_loc, a_loc] = gen_rand_state(Evt, names, opt, seed)
    arguments
        Evt
        names struct
        opt struct
        seed (1, 1) {mustBeNonnegative} = 0;
    end
    if seed~=0; rng(seed); end

    % Dimensions
    n = length(names.v);
    m = length(names.a);
    k = length(names.t);

    % Set location of each agent
    a_loc = cell(1, m);
    for i = 1:m
        a_loc{i} = names.v(randi([1 n]));
    end

    % Assign pending tasks
    t_loc = cell(1, k);
    for i = 1:k

        % Get location numbers where task i can be assigned
        Evti = Evt(:, i);
        idx = find(Evti == 1); % Location numbers where task i can be assigned
        num_t = length(idx); % Number of locations where the task can be assigned
        num_assigned = max(min(randi([0 num_t]), opt.max_tv), opt.min_tv); % Assign a random number of tasks
        cell_t = cell(1, num_assigned);

        % Assign tasks to random locations (allways in locations where task i can be assigned)
        if num_assigned > 0
            idx_rand = idx(randperm(num_t));
            idx_assigned = idx_rand(1:num_assigned);
            for j = 1:num_assigned
                cell_t(j) = names.v(idx_assigned(j));
            end
        end
        t_loc{i} = cell_t;
    end

end