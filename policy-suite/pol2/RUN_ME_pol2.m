% RUN_ME_POL2  Entry point for pol2, Policy 2 in the manuscript.
%
% This folder carries its own copy of calibration_result_<SPEC>.mat.
%
% Run from this folder, or with it on the path. The steps, in order:
%
%   setup_path()                   adds engine and pols
%   locate_corner_pol2(SPEC)       the regime-flip threshold in closed
%                                  form, plus a grid over the interior
%                                  tradeoff past it. 
%   welfare_search_pol(SPEC)       the real search, and the source of the
%                                  answer here rather than a cross-check
%   compute_policy_tables_pol2     the two comparison tables
%   print_pol_run_summary          the closing panel
%
% SPEC is set below; baseline by default.

clearvars
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setup_path.m'));

SPEC = 'baseline';

fprintf('=== pol2 [%s] (tau_K only; tau_R=0, tau_L/w_min FIXED) -- SPEC = %s ===\n\n', ...
        manuscript_policy_map('pol2'), SPEC);

% Step 1: the closed-form flip threshold and a grid search, which does
% not by itself establish an optimum.
locate_corner_pol2(SPEC);

% Step 2: the answer proper. This needs the Global Optimization Toolbox,
% so it will not run under Octave.
welfare_search_pol('pol2', SPEC);

% Step 3: the tables, preferring the search result to the grid candidate.
[T1, T2, rep0, rep1] = compute_policy_tables_pol2(SPEC);

% Step 4: the closing panel. No wage floor check, the floor being fixed
% throughout this search.
print_pol_run_summary('pol2', SPEC, rep0, rep1);
