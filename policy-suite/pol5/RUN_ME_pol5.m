% RUN_ME_POL5  Entry point for pol5, Policy 5 in the manuscript.
%
% This folder carries its own copy of calibration_result_<SPEC>.mat.
%
% Run from this folder, or with it on the path. The steps, in order:
%
%   setup_path()                       adds engine and pols
%   locate_corner_pol5(SPEC)           the threshold wage and the flat
%                                      interval below it, in closed form.
%                                      pol5 has no welfare search: its
%                                      optimum is analytic, so there is
%                                      no annealing cross-check to run
%   check_wmin_identification_pol      whether the reported floor is a real
%                                      optimum
%   compute_policy_tables_pol5         the two comparison tables
%   print_pol_run_summary              the closing panel
%
%
% SPEC is set below; baseline by default.

clearvars
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setup_path.m'));

SPEC = 'baseline';

fprintf('=== pol5 [%s] (w_min only; tau_R=0, tau_K/tau_L FIXED) -- SPEC = %s ===\n\n', ...
        manuscript_policy_map('pol5'), SPEC);

% Step 1: the threshold and the plateau below it, both exact.
locate_corner_pol5(SPEC);

% Step 2: the identification check. 
wmin_result = check_wmin_identification_pol('pol5', SPEC, 0);

% Step 3: the tables.
[T1, T2, rep0, rep1] = compute_policy_tables_pol5(SPEC, [], '', wmin_result);

% Step 4: the closing panel.
print_pol_run_summary('pol5', SPEC, rep0, rep1, wmin_result);
