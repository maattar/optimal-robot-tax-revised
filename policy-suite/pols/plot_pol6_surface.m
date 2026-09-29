function [Zc, vR_plot, vW_plot, marks] = plot_pol6_surface(SPEC, N_PLOT, opts)
% PLOT_POL6_SURFACE  Figure 5 in the paper
%
%   plot_pol6_surface
%   plot_pol6_surface('baseline', 25)
%   plot_pol6_surface('baseline', 25, opts)
%   [Zc, vR_plot, vW_plot, marks] = plot_pol6_surface(...)
if nargin < 1 || isempty(SPEC),   SPEC = 'baseline'; end
if nargin < 2 || isempty(N_PLOT), N_PLOT = 25;       end
if nargin < 3,                    opts = struct();   end

DEF.style             = 'surf';
DEF.tau_R_lim         = [];
DEF.w_min_lim         = [];
DEF.view              = [-55.56 30];   % the angle the manuscript figure uses
DEF.shading           = 'faceted';
DEF.zlim              = [];
DEF.face_color        = [0.92 0.92 0.92];
DEF.edge_color        = [0 0 0];
DEF.edge_width        = 1;
DEF.markers           = true;
DEF.pol1_marker       = true;
DEF.pol1_dir          = '';
DEF.show_legend       = true;
DEF.allow_sa_fallback = false;
DEF.save_stem         = '';
DEF.save_png          = '';
OPTIONAL_EMPTY = {'tau_R_lim', 'w_min_lim', 'zlim', 'save_png', 'save_stem', 'pol1_dir'};
fn = fieldnames(DEF);
for k = 1:numel(fn)
    if ~isfield(opts, fn{k}) || isempty(opts.(fn{k}))
        if ~any(strcmp(fn{k}, OPTIONAL_EMPTY)) || ~isfield(opts, fn{k})
            opts.(fn{k}) = DEF.(fn{k});
        end
    end
end
if isscalar(N_PLOT), N_PLOT = [N_PLOT N_PLOT]; end

    
if isempty(opts.pol1_dir)
    opts.pol1_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'pol1');
end

cfg = pol_config('pol6');
if ~isequal(cfg.free, {'tau_R', 'w_min'})
    error(['plot_pol6_surface: pol_config(''pol6'').free = {%s}, expected ' ...
           '{tau_R, w_min}. This figure hardcodes tau_R on x and w_min on ' ...
           'y. If the definition of pol6 has changed, this file has to ' ...
           'change with it, or an axis will be mislabelled.'], strjoin(cfg.free, ', '));
end

[~, res_file, corner_file] = pol_filenames(cfg, SPEC);

% LOAD 
if ~exist(res_file, 'file')
    error(['plot_pol6_surface: %s not found. Run welfare_search_pol(''pol6'', ' ...
           '''%s'') first -- this script draws that search''s Phase-1 grid ' ...
           'and never solves anything itself.'], res_file, SPEC);
end
S = load(res_file);

req = {'X_flat', 'SWF_flat', 'SWF_star', 'SWF_bench', 'cp'};
for k = 1:numel(req)
    if ~isfield(S, req{k})
        error(['plot_pol6_surface: %s lacks field ''%s''.'], res_file, req{k});
    end
end

i_tR = find(strcmp(cfg.free, 'tau_R'));
i_wm = find(strcmp(cfg.free, 'w_min'));

v_tR = unique(S.X_flat(:, i_tR), 'sorted');
v_wm = unique(S.X_flat(:, i_wm), 'sorted');
if numel(v_tR) * numel(v_wm) ~= size(S.X_flat, 1)
    error(['plot_pol6_surface: grid reconstruction failed -- %d x %d does not ' ...
           'match the %d saved points.'], numel(v_tR), numel(v_wm), size(S.X_flat, 1));
end

SWF2 = reshape(S.SWF_flat, [numel(unique(S.X_flat(:,1))), numel(unique(S.X_flat(:,2)))]);
if i_tR ~= 1
    SWF2 = SWF2.';
end
SWF2(~isfinite(SWF2)) = NaN;

n_fail = sum(~isfinite(S.SWF_flat));

if isfield(S, 'baseline_w_min')
    w_min_calib = S.baseline_w_min;
elseif isfield(S.cp, 'w_min')
    w_min_calib = S.cp.w_min;
else
    w_min_calib = NaN;
end

fprintf('\n===========================================================\n');
fprintf('  PLOT_POL6_SURFACE  SPEC=''%s''\n', SPEC);
fprintf('===========================================================\n');
fprintf('  Computation grid (from %s): %d x %d = %d points', ...
        res_file, numel(v_tR), numel(v_wm), numel(S.SWF_flat));
if n_fail > 0
    fprintf('  (%d did not converge -> holes in the surface)', n_fail);
end
fprintf('\n');
fprintf('  tau_R in [%.4f, %.4f]   w_min in [%.4f, %.4f]\n', ...
        v_tR(1), v_tR(end), v_wm(1), v_wm(end));


keepR = true(size(v_tR));
keepW = true(size(v_wm));
if ~isempty(opts.tau_R_lim)
    keepR = v_tR >= opts.tau_R_lim(1) & v_tR <= opts.tau_R_lim(2);
end
if ~isempty(opts.w_min_lim)
    keepW = v_wm >= opts.w_min_lim(1) & v_wm <= opts.w_min_lim(2);
end
if sum(keepR) < 2 || sum(keepW) < 2
    error(['plot_pol6_surface: Widen tau_R_lim / w_min_lim.'], ...
          sum(keepR), sum(keepW));
end
v_tR = v_tR(keepR);
v_wm = v_wm(keepW);
SWF2 = SWF2(keepR, keepW);

[idx_R, vR_plot] = local_decimate(v_tR, N_PLOT(1));
[idx_W, vW_plot] = local_decimate(v_wm, N_PLOT(2));
Zc = SWF2(idx_R, idx_W);

fprintf('  Plotting grid: %d x %d  (stride %d x %d; every plotted height is a ', ...
        numel(vR_plot), numel(vW_plot), ...
        max(1, round(numel(v_tR)/numel(vR_plot))), ...
        max(1, round(numel(v_wm)/numel(vW_plot))));
fprintf('real solved point)\n');
n_hole = sum(~isfinite(Zc(:)));
if n_hole > 0
    fprintf('  %d / %d plotted nodes are NaN (non-converged) -- gaps in the mesh.\n', ...
            n_hole, numel(Zc));
end


mk.sq.tau_R = 0;
mk.sq.w_min = w_min_calib;
mk.sq.SWF   = S.SWF_bench;
mk.sq.src   = sprintf('benchmark from %s', res_file);
mk.sq.ok    = isfinite(w_min_calib);


mk.p6.tau_R = NaN;
mk.p6.w_min = NaN;
mk.p6.SWF   = NaN;
mk.p6.src   = '';
mk.p6.ok    = true;
if exist(corner_file, 'file')
    C6 = load(corner_file);
    need6 = {'corner', 'SWF_corner'};
    miss6 = need6(~isfield(C6, need6));
    if ~isempty(miss6)
        error(['plot_pol6_surface: %s lacks %s. Re-run locate_corner_pol(''pol6'', ''%s'').'], ...
              corner_file, strjoin(miss6, ' and '), SPEC);
    end
    mk.p6.tau_R = C6.corner;
    mk.p6.w_min = 0;           
    mk.p6.SWF   = C6.SWF_corner;
    mk.p6.src   = sprintf('exact corner from %s', corner_file);
elseif opts.allow_sa_fallback
    sa_req = {[cfg.free{i_tR} '_star'], [cfg.free{i_wm} '_star']};
    if ~all(isfield(S, sa_req))
        error(['plot_pol6_surface: %s has neither a corner file nor the ' ...
               'annealing fields %s, so Policy 6''s marker cannot be placed ' ...
               'at all.'], res_file, strjoin(sa_req, ' and '));
    end
    mk.p6.tau_R = S.(sa_req{1});
    mk.p6.w_min = S.(sa_req{2});
    mk.p6.SWF   = S.SWF_star;
    mk.p6.src   = 'SA (grid-quantised)';
    warning(['plot_pol6_surface: %s not found, so Policy 6''s marker is the ' ...
             'annealing point. Welfare is flat along both of pol6''s ' ...
             'instruments at the optimum, so that point is an arbitrary ' ...
             'place on the plateau. Exploratory only -- do not use this ' ...
             'figure in the manuscript.'], corner_file);
else
    error(['plot_pol6_surface: %s not found. Run locate_corner_pol(''pol6'', ' ...
           '''%s'') first. Falling back on the annealing point would put ' ...
           'Policy 6''s marker at an arbitrary point of the welfare ' ...
           'plateau; pass opts.allow_sa_fallback = true if that is ' ...
           'genuinely what you want for a throwaway plot.'], corner_file, SPEC);
end

mk.p1.ok = false;
mk.p1.src = '(not drawn)';
if opts.pol1_marker
    cfg1 = pol_config('pol1');
    [~, ~, corner_file1] = pol_filenames(cfg1, SPEC);
    pol1_path = fullfile(opts.pol1_dir, corner_file1);
    if ~exist(pol1_path, 'file')
        error(['plot_pol6_surface: %s not found. Policy 1''s optimum is now ' ...
               'one of this figure''s three markers, so pol1 has to have ' ...
               'been run for the same SPEC before pol6''s figure can be ' ...
               'drawn: run RUN_ME_pol1 (or locate_corner_pol1(''%s'') from ' ...
               'the pol1 folder). Set opts.pol1_dir if pol1''s results live ' ...
               'somewhere else, or opts.pol1_marker = false for the old ' ...
               'two-marker figure.'], pol1_path, SPEC);
    end
    C1 = load(pol1_path);
    need1 = {'corner', 'SWF_corner', 'w_min_calib', 'cp'};
    miss1 = need1(~isfield(C1, need1));
    if ~isempty(miss1)
        error(['plot_pol6_surface: %s lacks %s. Re-run locate_corner_pol1 ' ...
               'with the current code.'], pol1_path, strjoin(miss1, ', '));
    end
    if isfield(C1, 'SPEC') && ~strcmp(C1.SPEC, SPEC)
        error(['plot_pol6_surface: SPEC MISMATCH. This figure is being drawn ' ...
               'for ''%s'' but %s identifies itself as ''%s''. Refusing to ' ...
               'put an optimum from one calibration on another''s ' ...
               'surface.'], SPEC, pol1_path, C1.SPEC);
    end

    CP_TOL = 1e-12;
    shared = {'alpha', 'psi', 'lam_bar', 'eta_bar', 'delta', 'theta', 'nu', ...
              'sigma', 'g_A', 'g_N', 'rho', 's_a1', 's_a2', 's_n', ...
              'r', 'tau_K', 'tau_L'};
    [gap_max, gap_name] = local_cp_gap(C1.cp, S.cp, shared);
    if ~isfinite(gap_max)
        error(['plot_pol6_surface: cannot compare pol1''s and pol6''s ' ...
               'parameters -- field ''%s'' is missing from one of them. ' ...
               'Both corner and search results have to come from the same ' ...
               'suite version.'], gap_name);
    end
    if gap_max > CP_TOL
        error(['plot_pol6_surface: pol1 and pol6 do not share a surface. ' ...
               'Field ''%s'' differs by %.3e (tolerance %.1e). pol1''s ' ...
               'optimum is a point of pol6''s welfare surface only while ' ...
               'the two hold the same parameters and the same fixed ' ...
               'instruments; they do not here, so the marker would be ' ...
               'misplaced. Re-run both from the same calibration.'], ...
              gap_name, gap_max, CP_TOL);
    end

    if isfinite(w_min_calib) && abs(C1.w_min_calib - w_min_calib) > CP_TOL
        error(['plot_pol6_surface: pol1''s calibrated floor (%.12f) is not ' ...
               'pol6''s (%.12f). The two markers would sit on different ' ...
               'w_min lines and the figure would not say what it claims ' ...
               'to.'], C1.w_min_calib, w_min_calib);
    end

    mk.p1.tau_R = C1.corner;
    mk.p1.w_min = C1.w_min_calib;
    mk.p1.SWF   = C1.SWF_corner;
    mk.p1.src   = sprintf('exact corner from %s (cp agrees with pol6 to %.1e)', ...
                          pol1_path, max(gap_max, 0));
    mk.p1.ok    = true;
end

fprintf('\n  --- Markers (all from solved optima, never off the grid) ---\n');
fprintf('  Status quo   : tau_R =%.9f  w_min =%.9f  SWF =%.10f\n', ...
        mk.sq.tau_R, mk.sq.w_min, mk.sq.SWF);
fprintf('                 [%s]\n', mk.sq.src);
if mk.p1.ok
    fprintf('  %s opt : tau_R*=%.9f  w_min =%.9f  SWF*=%.10f\n', ...
            manuscript_policy_map('pol1'), mk.p1.tau_R, mk.p1.w_min, mk.p1.SWF);
    fprintf('                 [%s]\n', mk.p1.src);
end
fprintf('  %s opt : tau_R*=%.9f  w_min*=%.9f  SWF*=%.10f\n', ...
        manuscript_policy_map('pol6'), mk.p6.tau_R, mk.p6.w_min, mk.p6.SWF);
fprintf('                 [%s]\n', mk.p6.src);


Z_RANGE = max(SWF2(:)) - min(SWF2(:));
local_surface_check('Status quo',                            mk.sq, v_tR, v_wm, SWF2, Z_RANGE);
if mk.p1.ok
    local_surface_check([manuscript_policy_map('pol1') ' optimum'], mk.p1, v_tR, v_wm, SWF2, Z_RANGE);
end
local_surface_check([manuscript_policy_map('pol6') ' optimum'],     mk.p6, v_tR, v_wm, SWF2, Z_RANGE);
fprintf('\n');

% =========================================================================
% DRAW
% =========================================================================
figure('Name', sprintf('pol6: SWF over (tau_R, w_min) -- SPEC=%s', SPEC), ...
       'Color', 'w', 'Position', [80 80 980 760]);

[Wg, Rg] = meshgrid(vW_plot, vR_plot);   % x = tau_R (Rg), y = w_min (Wg)

switch lower(opts.style)
    case 'meshc', hs = meshc(Rg, Wg, Zc);
    case 'mesh',  hs = mesh(Rg, Wg, Zc);
    case 'surf',  hs = surf(Rg, Wg, Zc);
    case 'surfc', hs = surfc(Rg, Wg, Zc);
    case 'surfl', hs = surfl(Rg, Wg, Zc);
    otherwise
        error(['plot_pol6_surface: unknown style ''%s''. Valid: meshc, mesh, ' ...
               'surf, surfc, surfl.'], opts.style);
end
try
    shading(opts.shading);
catch
    
end


try
    set(hs(1), 'FaceColor', opts.face_color, ...
               'EdgeColor', opts.edge_color, ...
               'LineWidth', opts.edge_width);
catch
    
end

COL_P1 = [0.00 0.35 0.85];


hleg = {}; sleg = {};
hold on;
if opts.markers
    if mk.sq.ok
        h = plot3(mk.sq.tau_R, mk.sq.w_min, mk.sq.SWF, 'o', 'MarkerSize', 15, ...
                  'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
        hleg{end+1} = h; sleg{end+1} = 'Status quo';
    end
    if mk.p1.ok
        h = plot3(mk.p1.tau_R, mk.p1.w_min, mk.p1.SWF, 'o', 'MarkerSize', 15, ...
                  'MarkerFaceColor', COL_P1, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        hleg{end+1} = h; sleg{end+1} = [manuscript_policy_map('pol1') ' optimum'];
    end
    h = plot3(mk.p6.tau_R, mk.p6.w_min, mk.p6.SWF, 'o', 'MarkerSize', 15, ...
              'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    hleg{end+1} = h; sleg{end+1} = [manuscript_policy_map('pol6') ' optimum'];
end
hold off;

xlim([-0.1 1])
xlabel('Robot tax');
ylabel('Minimum wage');
zlabel('Social welfare','FontWeight','bold');
view(opts.view(1), opts.view(2));
if ~isempty(opts.zlim), zlim(opts.zlim); end
try
    set(gca, 'Projection', 'perspective');
catch
    
end
set(gca, 'Box', 'off');
grid on;


if opts.markers && opts.show_legend && ~isempty(hleg)
    try
       
        legend([hleg{:}], sleg, 'Location', 'northeast', 'FontSize', 10, ...
               'AutoUpdate', 'off');
    catch
     
    end
end

local_in_view(gca, 'Status quo', mk.sq);
if mk.p1.ok, local_in_view(gca, [manuscript_policy_map('pol1') ' optimum'], mk.p1); end
local_in_view(gca, [manuscript_policy_map('pol6') ' optimum'], mk.p6);

if ~isempty(opts.save_png)
    print(gcf, '-dpng', '-r200', opts.save_png);
    fprintf('  Wrote %s\n', opts.save_png);
end
if ~isempty(opts.save_stem)
    
    set(gcf, 'PaperPositionMode', 'auto');
    pp = get(gcf, 'PaperPosition');
    set(gcf, 'PaperSize', pp(3:4));
    try
        savefig(gcf, [opts.save_stem '.fig']);
        fprintf('  Wrote %s.fig\n', opts.save_stem);
    catch
        fprintf('  (savefig unavailable here; no .fig written)\n');
    end
    print(gcf, '-depsc2', '-r300', [opts.save_stem '.eps']);
    fprintf('  Wrote %s.eps\n', opts.save_stem);
    print(gcf, '-dpdf', '-r300', [opts.save_stem '.pdf']);
    fprintf('  Wrote %s.pdf\n', opts.save_stem);
end


marks = mk;
tau_R_star = mk.p6.tau_R;
w_min_star = mk.p6.w_min;
SWF_star   = mk.p6.SWF;
corner_src = mk.p6.src;
SWF_bench  = mk.sq.SWF;
if mk.p1.ok
    tau_R_pol1 = mk.p1.tau_R;
    w_min_pol1 = mk.p1.w_min;
    SWF_pol1   = mk.p1.SWF;
else
    tau_R_pol1 = NaN;  w_min_pol1 = NaN;  SWF_pol1 = NaN;
end
pol1_src = mk.p1.src;

outfile = sprintf('pol6_surface_%s.mat', SPEC);
save(outfile, 'SPEC', 'Zc', 'vR_plot', 'vW_plot', 'marks', ...
     'tau_R_star', 'w_min_star', 'SWF_star', 'corner_src', ...
     'tau_R_pol1', 'w_min_pol1', 'SWF_pol1', 'pol1_src', ...
     'SWF_bench', 'w_min_calib', 'opts', 'N_PLOT');
fprintf('Saved %s\n\n', outfile);

end  




function [idx, v_out] = local_decimate(v, n_target)
    n = numel(v);
    if n_target >= n
        idx = (1:n).';
    else
        stride = max(1, round(n / n_target));
        idx = unique([1:stride:n, n]).';
    end
    v_out = v(idx);
end



function [gap_max, gap_name] = local_cp_gap(cpA, cpB, names)
    gap_max = 0; gap_name = '';
    for i = 1:numel(names)
        f = names{i};
        if ~isfield(cpA, f) || ~isfield(cpB, f)
            gap_max = Inf; gap_name = f;
            return
        end
        a = cpA.(f); b = cpB.(f);
        scale = max(1, max(abs(a), abs(b)));
        g = abs(a - b) / scale;
        if g > gap_max
            gap_max = g; gap_name = f;
        end
    end
    if isempty(gap_name), gap_name = names{1}; end
end



function local_surface_check(name, m, v_tR, v_wm, SWF2, z_range)
    if ~isfinite(m.tau_R) || ~isfinite(m.w_min) || ~isfinite(m.SWF)
        return
    end
    inside = m.tau_R >= v_tR(1) && m.tau_R <= v_tR(end) && ...
             m.w_min >= v_wm(1) && m.w_min <= v_wm(end);
    if ~inside
        fprintf('  %-18s : outside the computed grid, no surface comparison.\n', name);
        return
    end
    z_hat = NaN;
    try
        z_hat = interp2(v_wm(:).', v_tR(:), SWF2, m.w_min, m.tau_R, 'linear');
    catch
    end
    if ~isfinite(z_hat)
        [~, iR] = min(abs(v_tR - m.tau_R));
        [~, iW] = min(abs(v_wm - m.w_min));
        z_hat = SWF2(iR, iW);
        tagged = 'nearest node';
    else
        tagged = 'bilinear';
    end
    if ~isfinite(z_hat)
        fprintf('  %-18s : surface is NaN there, no comparison.\n', name);
        return
    end
    gap = abs(m.SWF - z_hat);
    rel = gap / max(z_range, eps);
    if rel > 0.02
        fprintf(['  %-18s : surface says %.10f, marker says %.10f (gap %.3e, ' ...
                 '%.1f%% of the z-range, %s) -- LARGE, worth a look.\n'], ...
                name, z_hat, m.SWF, gap, 100*rel, tagged);
    else
        fprintf('  %-18s : agrees with the surface to %.3e (%s).\n', name, gap, tagged);
    end
end




function local_in_view(ax, name, m)
    if ~isfinite(m.tau_R) || ~isfinite(m.w_min), return, end
    xl = get(ax, 'XLim'); yl = get(ax, 'YLim');
    if m.tau_R < xl(1) || m.tau_R > xl(2) || m.w_min < yl(1) || m.w_min > yl(2)
        warning(['plot_pol6_surface: the %s marker at (%.6f, %.6f) falls ' ...
                 'outside the drawn axes (x [%.3f, %.3f], y [%.3f, %.3f]) ' ...
                 'and will not appear. Widen the limits before using this ' ...
                 'figure.'], name, m.tau_R, m.w_min, xl(1), xl(2), yl(1), yl(2));
    end
end
