% RUN_ME_POL1  Entry point for pol1, Policy 1 in the manuscript.
%
% This folder carries its own copy of calibration_result_<SPEC>.mat, taken
% from the shared calibration used throughout the suite.
%
% Run from this folder, or with it on the path. The steps, in order:
%
%   setup_path()                   adds engine and pols
%   locate_corner_pol1(SPEC)       the exact corner, in closed form
%   welfare_search_pol(SPEC)       a grid pass then simulated annealing,
%                                  here a cross-check on the corner rather
%                                  than the source of the answer
%   compute_policy_tables_pol1     the two comparison tables
%   print_pol_run_summary          the closing panel
%   plot_pol1_sweep                the manuscript figure, six panels
%                                  against the robot tax, written to
%                                  figure_tauR.fig, .eps and .pdf in this
%                                  folder
%
% SPEC is set below. Options: baseline, auto33, equal_shares,
% acemoglu2020. 

clearvars
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setup_path.m'));

SPEC    = 'baseline';
N_SWEEP = 400;          % points in the figure's robot-tax sweep

fprintf('=== pol1 [%s] (tau_R only; tau_K/tau_L/w_min FIXED) -- SPEC = %s ===\n\n', ...
        manuscript_policy_map('pol1'), SPEC);

% Step 1: the closed-form corner, which is the optimum for pol1, welfare
% being flat above it on the same argument that settles the pol6 corner.
locate_corner_pol1(SPEC);

% Step 2: an optional cross-check by simulated annealing. This needs the
% Global Optimization Toolbox, so it will not run under Octave. It can be
% skipped where only the closed-form result is wanted.
welfare_search_pol('pol1', SPEC);

% Step 3: the tables, preferring the exact corner to the search result.
[T1, T2, rep0, rep1] = compute_policy_tables_pol1(SPEC);

% Step 4: the closing panel. No wage floor check, the floor being fixed
% throughout this search as it is under pol4.
print_pol_run_summary('pol1', SPEC, rep0, rep1);

% Step 5: the manuscript figure. The stem carries no SPEC -- drawing an 
% alternative specification here overwrites the baseline figure.
plot_pol1_sweep(SPEC, N_SWEEP, [], 'figure_tauR');
