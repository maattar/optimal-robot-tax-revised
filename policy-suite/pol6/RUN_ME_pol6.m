% RUN_ME_POL6  Entry point for pol6, Policy 6 in the manuscript.
%
% This folder carries its own copy of calibration_result_<SPEC>.mat.
%
% RUN POL1 FIRST. The last step draws the manuscript figure, which carries
% a marker at Policy 1's optimum and reads it from pol1's own corner
% result, so pol1 must have been run for the same SPEC. Running out of
% order stops at that step with a message saying so.
%
% Run from this folder, or with it on the path. The steps, in order:
%
%   setup_path()                       adds engine and pols
%   welfare_search_pol(SPEC)           a grid pass then simulated
%                                      annealing, which for pol6 is a
%                                      check on the plateau rather than
%                                      the answer
%   locate_corner_pol(SPEC)            the exact corner, in closed form
%   check_wmin_identification_pol      whether the reported floor is a real
%                                      optimum.
%   plot_policy_grid_pol               the two tables, and a working
%                                      surface plot that is not the
%                                      manuscript figure
%   print_pol_run_summary              the closing panel
%   plot_pol6_surface                  the manuscript figure, three
%                                      markers, written to wmintaur.fig,
%                                      .eps and .pdf in this folder
%
% SPEC is set below; baseline by default.

clearvars
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setup_path.m'));

SPEC = 'baseline';

fprintf('=== pol6 [%s] (tau_R, w_min) -- SPEC = %s ===\n\n', ...
        manuscript_policy_map('pol6'), SPEC);

welfare_search_pol('pol6', SPEC);
locate_corner_pol('pol6', SPEC);
wmin_result = check_wmin_identification_pol('pol6', SPEC);
[T1, T2, rep0, rep1] = plot_policy_grid_pol('pol6', SPEC, wmin_result);
print_pol_run_summary('pol6', SPEC, rep0, rep1, wmin_result);

% N_PLOT is the drawn mesh only; welfare is read off the search's own
% 636 x 636 grid.
fig_opts.save_stem = '3dfigfinal';
plot_pol6_surface(SPEC, 50, fig_opts);
