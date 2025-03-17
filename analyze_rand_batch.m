%% This script is used to analyze the results obtained from run_rand_grid()
clear; clc;

%% Determine which results to analyze
get_n = 72;

%% Load files
mat = dir('results/*.mat');

%% Determine number of elegible tests
num_tests = 0;
load_names = cell(1, 1);
for i = 1:length(mat)

    % Load test
    test = load(['results/' mat(i).name]);
    
    if isnan(get_n) || (test.test_info.n == get_n)
        num_tests = num_tests + 1;
        load_names{num_tests} = mat(i).name;
    end

end

%% Initialize matrices containing the results of eligible tests
if num_tests == 0
    error("No tests loaded. No tests in results/ satisfy the criteria. Aborting")
end
test = load(['results/' load_names{1}]);
fn = fieldnames(test.res);
res = struct();
for j = 1:length(fn)
    res.(fn{j}) = [];
end
% res = rmfield(res, 'field_name'); % To remove a field
fn = fieldnames(res);

%% Load data
for i = 1:num_tests

    % Load test
    test = load(['results/' load_names{i}]);

    % Save results into res
    for j = 1:length(fn)
        if isfield(test.res, fn{j})
            if ~isnan(test.res.(fn{j}))
                res.(fn{j})(end+1) = test.res.(fn{j});
            end
        end
    end

end

%% Analyze results
time_avrg_c = mean(res.time_c);
time_med_c = median(res.time_c);
time_max_c = max(res.time_c);
time_min_c = min(res.time_c);
time_avrg_d = mean(res.time_d_par_total);
time_med_d = median(res.time_d_par_total);
time_max_d = max(res.time_d_par_total);
time_min_d = min(res.time_d_par_total);
err_V_avrg = mean(res.err_V);
err_V_med = median(res.err_V);
err_V_max = max(res.err_V);
err_V_min = min(res.err_V);
k_d_avrg = mean(res.k_d);
k_d_med = median(res.k_d);
k_d_max = max(res.k_d);
k_d_min = min(res.k_d);
eflag_1_d = sum(res.e_flag_d == 1);
eflag_2_d = sum(res.e_flag_d == 2);
eflag_m1_d = sum(res.e_flag_d == -1);
eflag_1_c = sum(res.e_flag_c == 1);

%% Print results
fprintf("*** REPORT ***\n")
fprintf("> Number of loaded tests: %d (n: %d)\n", num_tests, get_n);
fprintf("> Results ('C' for centralized, 'H' for heuristic)\n");
fprintf("   > Time Avrg.  \t> C: %1.6f \t H: %1.6f\n", time_avrg_c, time_avrg_d);
fprintf("   > Time Med.   \t> C: %1.6f \t H: %1.6f\n", time_med_c, time_med_d);
fprintf("   > Time Max.   \t> C: %1.6f \t H: %1.6f\n", time_max_c, time_max_d);
fprintf("   > Time Min.   \t> C: %1.6f \t H: %1.6f\n", time_min_c, time_min_d);
fprintf("   > Iter Avrg.  \t>          \t H: %1.6f\n", k_d_avrg);
fprintf("   > Iter Med.   \t>          \t H: %1.2f\n", k_d_med);
fprintf("   > Iter Max.   \t>          \t H: %d\n", k_d_max);
fprintf("   > Iter Min.   \t>          \t H: %d\n", k_d_min);
fprintf("   > Err. Avrg.  \t> %g\n", err_V_avrg);
fprintf("   > Err. Med.   \t> %g\n", err_V_med);
fprintf("   > Err. Max.   \t> %g\n", err_V_max);
fprintf("   > Err. Min.   \t> %g\n", err_V_min);
fprintf("   > Finihed     \t> C: %d. \t Flags H -> 1: %d, 2: %d, -1: %d\n", eflag_1_c, eflag_1_d, eflag_2_d, eflag_m1_d);