% RUN_ME_POL3  Entry point for pol3, Policy 3 in the manuscript.
%
% This folder carries its own copy of calibration_result_<SPEC>.mat.
%
% Run from this folder, or with it on the path. The steps, in order:
%
%   setup_path()                   adds engine and pols
%   locate_corner_pol3(SPEC)       pol2's flip threshold, which carries
%                                  over unchanged, plus a grid over the
%                                  surface in both rates
%   welfare_search_pol(SPEC)       the real search: tau_L added
%   compute_policy_tables_pol3     the two comparison tables
%   print_pol_run_summary          the closing panel
%
% SPEC is set below; baseline by default.

clearvars
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setup_path.m'));

SPEC = 'baseline';

fprintf('=== pol3 [%s] (tau_K, tau_L; tau_R=0, w_min FIXED at status quo) -- SPEC = %s ===\n\n', ...
        manuscript_policy_map('pol3'), SPEC);

% Step 1: the closed-form flip threshold, the same as pol2, and a search
% over a two-dimensional grid, which does not by itself establish an
% optimum.
locate_corner_pol3(SPEC);

% Step 2: the answer proper. This needs the Global Optimization Toolbox,
% so it will not run under Octave.
welfare_search_pol('pol3', SPEC);

% Step 3: the tables, preferring the search result to the grid candidate.
[T1, T2, rep0, rep1] = compute_policy_tables_pol3(SPEC);

% Step 4: the closing panel. No wage floor check, the floor being fixed
% throughout this search.
print_pol_run_summary('pol3', SPEC, rep0, rep1);
