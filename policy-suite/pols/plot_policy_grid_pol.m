function [T1, T2, rep0, rep1] = plot_policy_grid_pol(pol_id, SPEC, wmin_result)
% PLOT_POLICY_GRID_POL  Post-processing for pol6 and pol4
%
%   plot_policy_grid_pol('pol6')
%   plot_policy_grid_pol('pol4', 'baseline')
%   plot_policy_grid_pol('pol6', 'baseline', wmin_result)
%
% The manuscript figure is plot_pol6_surface, not the one drawn here.
%

if nargin < 2 || isempty(SPEC), SPEC = 'baseline'; end
if nargin < 3, wmin_result = []; end

cfg = pol_config(pol_id);
if ~cfg.built
    error('plot_policy_grid_pol: %s is not yet built.', pol_id);
end
nfree = numel(cfg.free);
if nfree ~= 2 && nfree ~= 3
    error(['plot_policy_grid_pol: %s has %d free instruments. Only pol6 ' ...
           '(2 free, gets a figure) and pol4 (3 free, tables only, no ' ...
           'figure per Part B.1) are wired up here.'], pol_id, nfree);
end
if ~any(strcmp(pol_id, {'pol6', 'pol4'}))
    error(['plot_policy_grid_pol: %s is not supported here. This file ' ...
           'hardcodes tau_R as cfg.free{1} throughout (see the "tau_R is ' ...
           'always cfg.free{1} so far" comment below), so a policy that ' ...
           'fixes tau_R would have its capital tax mislabelled as the ' ...
           'robot tax.'], pol_id);
end

[~, res_file, corner_file] = pol_filenames(cfg, SPEC);


% Load 
if ~exist(res_file, 'file')
    error(['plot_policy_grid_pol: %s not found. Run welfare_search_pol ' ...
           'first (same pol_id and SPEC).'], res_file);
end
S = load(res_file);
cp = S.cp;

x_star = zeros(1, nfree);
for i = 1:nfree
    x_star(i) = S.([cfg.free{i} '_star']);
end
SWF_star  = S.SWF_star;
SWF_bench = S.SWF_bench;

fprintf('\nLoaded %s.\n', res_file);
fprintf('  SA-reported: %s\n', join_star(cfg.free, x_star));
fprintf('  SWF* = %.8f    SWF(bench) = %.8f\n\n', SWF_star, SWF_bench);


if isfield(S, 'baseline_w_min'), w_min_calib = S.baseline_w_min; else, w_min_calib = cp.w_min; end
if isfield(S, 'baseline_tau_K'), tau_K_base = S.baseline_tau_K; elseif isfield(cp, 'tau_K'), tau_K_base = cp.tau_K; else, tau_K_base = NaN; end
if isfield(S, 'baseline_tau_L'), tau_L_base = S.baseline_tau_L; elseif isfield(cp, 'tau_L'), tau_L_base = cp.tau_L; else, tau_L_base = NaN; end


tau_R_star_sa = x_star(1);
tau_R_report  = tau_R_star_sa;
corner_src    = 'SA (grid-quantised -- run locate_corner_pol for the exact value)';
have_corner   = false;
if exist(corner_file, 'file')
    Cres = load(corner_file);
    if isfield(Cres, 'C_star') || isfield(Cres, 'corner')
        have_corner = true;
    end
end

if cfg.needs_euler
    tau_K_star = x_star(strcmp(cfg.free, 'tau_K'));
    if have_corner && isfield(Cres, 'C_star')
        C_star = Cres.C_star;
        fprintf('Loaded exact corner data from %s (C = %.12f).\n', corner_file, C_star);
    else
        fprintf('%s not found or lacks C_star; deriving C inline.\n', corner_file);
        r_probe = euler_closure(cp.Phi, tau_K_base);
        [~, C_star] = noauto_wage_pol(r_probe, cp);
        fprintf('  C = %.12f  (from w_a1^{c,NA} at tau_K_base)\n', C_star);
    end
end

if cfg.needs_euler
    r_at_star = euler_closure(cp.Phi, tau_K_star);
    corner_at_star = corner_of_pol(r_at_star, cp.delta, cp.lam_bar, w_min_calib, C_star, cp.alpha);
    if have_corner
        tau_R_report = corner_at_star;
        corner_src   = sprintf('exact corner at tau_K* (from %s)', corner_file);
    end
    fprintf('  tau_R^corner(tau_K*) = %.6f\n', corner_at_star);
    fprintf('  gap  tau_R*(SA) - corner = %+.6f  (plateau quantisation)\n', ...
            tau_R_star_sa - corner_at_star);
else
    if have_corner && isfield(Cres, 'corner')
        tau_R_report = Cres.corner;
        corner_src   = sprintf('exact corner (from %s)', corner_file);
    end
end
fprintf('  Tables below use tau_R = %.6f  [%s]\n\n', tau_R_report, corner_src);


% grid 
if ~isfield(S, 'X_flat')
    error(['plot_policy_grid_pol: %s has no X_flat field. Was it produced ' ...
           'by welfare_search_pol.m?'], res_file);
end
X_flat = S.X_flat;   % [N_TOTAL x nfree], as saved by welfare_search_pol.m
vecs = cell(1, nfree);
Ns   = zeros(1, nfree);
for i = 1:nfree
    vecs{i} = unique(S.X_flat(:, i), 'sorted');
    Ns(i)   = numel(vecs{i});
end
if prod(Ns) ~= size(S.X_flat, 1)
    error(['plot_policy_grid_pol: grid reconstruction failed -- prod(N) ' ...
           '(%d) does not match the %d saved points.'], prod(Ns), size(S.X_flat,1));
end
SWFnD = reshape(S.SWF_flat, Ns);
SWFnD(~isfinite(SWFnD)) = NaN;

n_fail = sum(~isfinite(S.SWF_flat));
if n_fail > 0
    fprintf('  Note: %d / %d grid points did not converge (blank on the contour panels).\n\n', ...
            n_fail, numel(S.SWF_flat));
end

% pol6 figure
if strcmp(pol_id, 'pol6')
    % Welfare over the robot tax and the wage floor.

    i_tR = find(strcmp(cfg.free, 'tau_R'));
    i_wm = find(strcmp(cfg.free, 'w_min'));
    v_tR = vecs{i_tR}; v_wm = vecs{i_wm};
    SWF2 = SWFnD;
    if i_tR ~= 1
        SWF2 = SWF2';  
    end

    
    N_PLOT_TARGET = 40;
    stride_R = max(1, round(numel(v_tR) / N_PLOT_TARGET));
    stride_W = max(1, round(numel(v_wm) / N_PLOT_TARGET));
    idx_R = [1:stride_R:numel(v_tR), numel(v_tR)]; idx_R = unique(idx_R);
    idx_W = [1:stride_W:numel(v_wm), numel(v_wm)]; idx_W = unique(idx_W);
    vR_plot = v_tR(idx_R); vW_plot = v_wm(idx_W);
    Z_plot  = SWF2(idx_R, idx_W);

    figure('Name', sprintf('Optimal %s policy -- SWF over (tau_R, w_min)', pol_id), ...
           'Color', 'w', 'Position', [100 100 900 700]);
    [Wg, Rg] = meshgrid(vW_plot, vR_plot);   % X=tau_R, Y=w_min convention below
    meshc(Rg, Wg, Z_plot);
    colormap('parula');
    hold on;
    x_bench = zeros(1, nfree);
    for i = 1:nfree
        if strcmp(cfg.free{i}, 'tau_R'), x_bench(i) = 0.0; else, x_bench(i) = w_min_calib; end
    end
 
    plot3(x_bench(i_tR), x_bench(i_wm), SWF_bench, 'o', 'MarkerSize', 10, ...
          'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', 'LineWidth', 1.0);
    plot3(tau_R_report, x_star(i_wm), SWF_star, 'o', 'MarkerSize', 10, ...
          'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth', 1.0);
    hold off;
    xlabel('\tau_R'); ylabel('w_{min}'); zlabel('SWF');
    view(3);
    set(gca, 'Projection', 'perspective', 'Box', 'off');
    grid on;
    title(sprintf(['SWF over (\\tau_R, w_{min}), plotted on a coarsened %d\\times%d grid ' ...
                   '(calc. grid: %d\\times%d). Black ball = status quo, red ball = optimum.'], ...
                  numel(vR_plot), numel(vW_plot), numel(v_tR), numel(v_wm)), ...
          'FontWeight', 'normal', 'FontSize', 10);

else
  
    fprintf(['  (No figure for %s -- figures are not economically legible ' ...
             'beyond 2 free instruments. Tables only, below.)\n\n'], pol_id);
end


PLATEAU_EPS = 1e-6;
x_bench = zeros(1, nfree);
x_opt   = x_star;
for i = 1:nfree
    switch cfg.free{i}
        case 'tau_R', x_bench(i) = 0.0;         x_opt(i) = tau_R_report + PLATEAU_EPS;
        case 'tau_K', x_bench(i) = tau_K_base;
        case 'tau_L', x_bench(i) = tau_L_base;
        case 'w_min', x_bench(i) = w_min_calib;
    end
end
rep0 = compute_policy_report_pol(x_bench, cfg, cp);
rep1 = compute_policy_report_pol(x_opt,   cfg, cp);
rep1.tau_R = tau_R_report;   % report the EXACT corner

if ~rep0.converged
    error('plot_policy_grid_pol: status-quo report failed to converge.');
end
if ~rep1.converged
    error('plot_policy_grid_pol: optimal report failed to converge.');
end


if any(strcmp('w_min', cfg.free)) && ~isempty(wmin_result) && wmin_result.is_indeterminate
    w_min_raw = rep1.w_min;
    if abs(w_min_raw) > 1e-9
        rep1.w_min = 0;
        fprintf(['  Note: w_min identification verdict is INDETERMINATE (flat plateau) -- ' ...
                 'reporting the recommended focal point w_min=0 in the tables below ' ...
                 '(raw search value was %.6f).\n\n'], w_min_raw);
    end
end


[T1, T2, cev] = print_policy_tables_pol(rep0, rep1, cfg, cp, pol_id, '');

% Save
outfile = sprintf('policy_grid_result_%s_%s.mat', cfg.tag, SPEC);
save(outfile, 'SPEC', 'vecs', 'SWFnD', 'rep0', 'rep1', 'T1', 'T2', 'cev', ...
     'x_star', 'tau_R_report', 'have_corner', 'w_min_calib', 'cp');
fprintf('\nSaved %s\n\n', outfile);

end  

function s = join_star(names, x)
    parts = cell(1, numel(names));
    for i = 1:numel(names)
        parts{i} = sprintf('%s*=%.6f', names{i}, x(i));
    end
    s = strjoin(parts, '  ');
end
