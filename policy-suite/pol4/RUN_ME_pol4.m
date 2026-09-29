% RUN_ME_POL4  Entry point for pol4, Policy 4 in the manuscript.
%
% This folder carries its own copy of calibration_result_<SPEC>.mat.
%
% Run from this folder, or with it on the path. The steps, in order:
%
%   setup_path()                   adds engine and pols
%   welfare_search_pol(SPEC)       runs first here, because the corner
%                                  locator reads tau_L* off its result
%   locate_corner_pol(SPEC)        the corner, a frontier in tau_K rather
%                                  than a single point
%   plot_policy_grid_pol           the two tables; no figure, three free
%                                  instruments
%   print_pol_run_summary          the closing panel
%
%
% SPEC is set below; baseline by default.

clearvars
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setup_path.m'));

SPEC = 'baseline';

fprintf('=== pol4 [%s] (tau_R, tau_K, tau_L) -- SPEC = %s ===\n\n', ...
        manuscript_policy_map('pol4'), SPEC);

welfare_search_pol('pol4', SPEC);
locate_corner_pol('pol4', SPEC);
[T1, T2, rep0, rep1] = plot_policy_grid_pol('pol4', SPEC);
print_pol_run_summary('pol4', SPEC, rep0, rep1);   
