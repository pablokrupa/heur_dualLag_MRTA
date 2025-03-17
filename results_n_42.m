%% Script used to reproduce the numerical results presented
clear; clc;
n_rows = 7; % Number of rows of the grid
n_cols = 6; % Number of columns of the grid
n = n_rows * n_cols;
m = n_rows;
k = n_cols;
seed.graph = n;
rng(seed.graph); % Use dimension of grid to determine RNG seed

%% Generate random seeds for assignment of locations and tasks (state)
num_tests = 10; % Number of times the grid is solved for different states
seed_state = randi([1, 1000], [1, num_tests]);
% This should be:
%   seed_state = [155,741,264,534,15,919,901,34,957,138];

%% Parameters of the grid
param = gen_param_grid(n_rows, n_cols, m, k);

%% Select other options
save_me = true;
save_name = 'grid';

%% Run tests
for i = 1:num_tests
    seed.state = seed_state(i);
    fprintf("*** Running test #%d of %d ***\n", i, num_tests);
    run_rand_grid(param, seed, save_me, save_name, save_add_time=false);
end
