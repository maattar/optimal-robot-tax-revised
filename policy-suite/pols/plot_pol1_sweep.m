function [G, rep0, rep1] = plot_pol1_sweep(SPEC, N_PTS, TAU_R_LIM, SAVE_STEM)
% PLOT_POL1_SWEEP  Figure 4
%
%   plot_pol1_sweep
%   plot_pol1_sweep('baseline', 400)
%   plot_pol1_sweep('baseline', 400, [-0.10 1.50])
%   plot_pol1_sweep('baseline', 400, [], 'figure_tauR')
%   [G, rep0, rep1] = plot_pol1_sweep(...)
%

if nargin < 1 || isempty(SPEC),      SPEC = 'baseline'; end
if nargin < 2 || isempty(N_PTS),     N_PTS = 250;       end
if nargin < 3,                       TAU_R_LIM = [];    end
if nargin < 4 || isempty(SAVE_STEM), SAVE_STEM = '';    end

    
PANELS = {
    'SWF',            'Social welfare',              1
    'u_a1',           'Unemployment',                1
    'avg_HtM',        'Avg. hand-to-mouth income',   1
    'redist_budget',  'Redistribution budget',       1
    'IRshare',        'Robotic investment share',    1
    'avg_Ricardian',  'Avg. Ricardian income',       1
};
NR = 2; NC = 3;
if size(PANELS, 1) ~= NR * NC
    error('plot_pol1_sweep: PANELS has %d rows but the layout is %dx%d.', ...
          size(PANELS, 1), NR, NC);
end

cfg = pol_config('pol1');
if ~isequal(cfg.free, {'tau_R'})
    error(['plot_pol1_sweep: pol_config(''pol1'').free = {%s}, expected ' ...
           '{tau_R} only. This figure hardcodes tau_R as the x-axis of ' ...
           'all six panels.'], ...
          strjoin(cfg.free, ', '));
end

[cal_file, res_file, corner_file] = pol_filenames(cfg, SPEC);


if ~exist(cal_file, 'file')
    error(['plot_pol1_sweep: %s not found. Run this from the pol1/ folder ' ...
           '(after copying the calibration .mat in from pol6/), or add ' ...
           'that folder to the path.'], cal_file);
end
S = load(cal_file);
if ~isfield(S, 'zeta_star') || ~isfield(S, 'fp')
    error('plot_pol1_sweep: %s must contain both ''zeta_star'' and ''fp''.', cal_file);
end
if isfield(S, 'SPEC') && ~isempty(SPEC) && ~strcmp(S.SPEC, SPEC)
    error(['plot_pol1_sweep: SPEC MISMATCH. Requested ''%s'' but %s ' ...
           'identifies itself as ''%s''.'], SPEC, cal_file, S.SPEC);
end
zeta_star = S.zeta_star;  fp = S.fp;

cp.alpha   = zeta_star(1);
cp.psi     = zeta_star(2);
cp.lam_bar = zeta_star(3);
cp.eta_bar = zeta_star(4);
w_min_calib = zeta_star(5);

cp.delta = fp.delta;  cp.theta = fp.theta;  cp.nu   = fp.nu;
cp.sigma = fp.sigma;  cp.g_A   = fp.g_A;    cp.g_N  = fp.g_N;
cp.rho   = fp.rho;    cp.s_a1  = fp.s_a1;   cp.s_a2 = fp.s_a2;
cp.s_n   = 1 - fp.s_a1 - fp.s_a2;
cp.verbose = false;

cp.r     = fp.r;         
cp.tau_K = fp.tau_K;
cp.tau_L = fp.tau_L;
cp.w_min = w_min_calib;   

fprintf('\n===========================================================\n');
fprintf('  PLOT_POL1_SWEEP  SPEC=''%s''  (tau_K, tau_L, w_min ALL FIXED)\n', SPEC);
fprintf('===========================================================\n');
fprintf('  alpha=%.6f  psi=%.6f  lam_bar=%.6f  eta_bar=%.6f\n', ...
        cp.alpha, cp.psi, cp.lam_bar, cp.eta_bar);
fprintf('  r (FIXED) = %.9f   w_min (FIXED) = %.9f\n', cp.r, w_min_calib);


tau_R_star = NaN;
corner_src = 'UNKNOWN (no corner_result_* and no welfare_result_* found)';
if exist(corner_file, 'file')
    Cres = load(corner_file);
    if isfield(Cres, 'corner')
        tau_R_star = Cres.corner;
        corner_src = sprintf('EXACT closed form, from %s', corner_file);
    end
end
if isnan(tau_R_star) && exist(res_file, 'file')
    R2 = load(res_file);
    if isfield(R2, 'tau_R_star')
        tau_R_star = R2.tau_R_star;
        corner_src = sprintf(['SA (grid-quantised) from %s -- run ' ...
                              'locate_corner_pol1 for the exact value'], res_file);
    end
end
fprintf('  tau_R* = %s\n', num2str(tau_R_star, '%.12f'));
fprintf('  source : %s\n', corner_src);


if isempty(TAU_R_LIM)
    if isnan(tau_R_star)
        TAU_R_LIM = [-0.10, 1.0];
    else
        TAU_R_LIM = [-0.10, max(1.0, 2 * tau_R_star)];
    end
end
tau_grid = linspace(TAU_R_LIM(1), TAU_R_LIM(2), N_PTS)';



PLATEAU_EPS = 1e-6;
if ~isnan(tau_R_star) && tau_R_star > TAU_R_LIM(1) && tau_R_star < TAU_R_LIM(2)
    tau_grid = unique([tau_grid; tau_R_star; tau_R_star + PLATEAU_EPS]);
end

if 0 > TAU_R_LIM(1) && 0 < TAU_R_LIM(2)
    tau_grid = unique([tau_grid; 0]);
end
N = numel(tau_grid);

fprintf('  sweep  : tau_R in [%.4f, %.4f], %d points\n\n', ...
        TAU_R_LIM(1), TAU_R_LIM(2), N);


nP = size(PANELS, 1);
Y  = nan(N, nP);
converged   = false(N, 1);
regime      = repmat({'?'}, N, 1);
task_regime = repmat({'?'}, N, 1);
auto_share  = nan(N, 1);

t0 = tic;
for i = 1:N
    try
        rep = compute_policy_report_pol(tau_grid(i), cfg, cp);
    catch ME
        regime{i} = sprintf('ERROR: %s', ME.message);
        continue
    end
    regime{i} = rep.regime;
    if ~rep.converged
        continue
    end
    converged(i)   = true;
    task_regime{i} = rep.task_regime;
    auto_share(i)  = rep.automation_share;
    for j = 1:nP
        Y(i, j) = rep.(PANELS{j, 1}) * PANELS{j, 3};
    end
end
elapsed = toc(t0);

n_ok = sum(converged);
fprintf('Sweep complete (%.1f s). Converged: %d / %d.\n', elapsed, n_ok, N);
if n_ok < N
    bad = find(~converged);
    fprintf('  NOT converged at %d point(s); tau_R = ', numel(bad));
    fprintf('%.4f ', tau_grid(bad(1:min(12, numel(bad)))));
    if numel(bad) > 12, fprintf('... (+%d more)', numel(bad) - 12); end
    fprintf('\n');
    fprintf(['  In Octave this is the known optimoptions/fsolve wall on the ' ...
             'automation branch, not a result. In MATLAB it is not expected ' ...
             '-- investigate before using the figure.\n']);
end
if n_ok == 0
    error(['plot_pol1_sweep: zero points converged -- nothing to draw. ' ...
           'This is the MATLAB-only solver wall if you are in Octave.']);
end


fprintf('  Regime census over converged points:\n');
uregs = unique(regime(converged));
for k = 1:numel(uregs)
    sel = converged & strcmp(regime, uregs{k});
    fprintf('    %-14s : %4d point(s), tau_R in [%.4f, %.4f]\n', ...
            uregs{k}, sum(sel), min(tau_grid(sel)), max(tau_grid(sel)));
end
utasks = unique(task_regime(converged));
fprintf('  Task-space sub-regimes present: %s\n\n', strjoin(utasks', ', '));



figure('Name', sprintf('pol1: policy sweep vs. tau_R (SPEC=%s)', SPEC), ...
       'Color', 'w', 'Position', [80 80 1250 720]);

for j = 1:nP
    subplot(NR, NC, j);
    plot(tau_grid, Y(:, j), '-', 'LineWidth', 1.6, 'Color', [0 0.35 0.70]);
    hold on;
    local_vline(0,          'k--', 1.1);            % status quo
    if ~isnan(tau_R_star)
        local_vline(tau_R_star, 'r-', 1.4);         % optimum
    end
    hold off;
    grid on;
    box off;
    xlim(TAU_R_LIM);
    % All six panels share the x-axis, so only the bottom row is labelled.
    if j > NC
        xlabel('Robot tax');
    end
    title(PANELS{j, 2}, 'FontWeight', 'bold');
end


if ~isempty(SAVE_STEM)
    
    set(gcf, 'PaperPositionMode', 'auto');
    pp = get(gcf, 'PaperPosition');
    set(gcf, 'PaperSize', pp(3:4));
    try
        savefig(gcf, [SAVE_STEM '.fig']);
        fprintf('  Wrote %s.fig\n', SAVE_STEM);
    catch
        fprintf('  (savefig unavailable here; no .fig written)\n');
    end
    print(gcf, '-depsc2', '-r300', [SAVE_STEM '.eps']);
    fprintf('  Wrote %s.eps\n', SAVE_STEM);
    print(gcf, '-dpdf', '-bestfit', '-r300', [SAVE_STEM '.pdf']);
    fprintf('  Wrote %s.pdf\n', SAVE_STEM);
end


rep0 = local_safe_report(0, cfg, cp);
if ~isnan(tau_R_star)
    rep1 = local_safe_report(tau_R_star + PLATEAU_EPS, cfg, cp);
    if ~isempty(rep1)
        rep1.tau_R = tau_R_star;   % report the EXACT corner
    end
else
    rep1 = [];
end

fprintf('  Panel values at the two marked points:\n');
fprintf('    %-26s %14s %14s\n', 'variable', 'tau_R = 0', 'tau_R = tau_R*');
for j = 1:nP
    v0 = local_pick(rep0, PANELS{j, 1}, PANELS{j, 3});
    v1 = local_pick(rep1, PANELS{j, 1}, PANELS{j, 3});
    fprintf('    %-26s %14.8f %14.8f\n', PANELS{j, 1}, v0, v1);
end
if isempty(rep0) || isempty(rep1)
    fprintf(['    (NaN column(s) above = that endpoint did not solve. In Octave ' ...
             'the tau_R=0 column is expected to be NaN: the status quo is in ' ...
             'the automation regime, which needs MATLAB.)\n']);
end
fprintf('\n');


G = struct();
G.SPEC        = SPEC;
G.tau_R       = tau_grid;
G.panels      = PANELS(:, 1)';
G.Y           = Y;
G.converged   = converged;
G.regime      = regime;
G.task_regime = task_regime;
G.automation_share = auto_share;
G.tau_R_star  = tau_R_star;
G.corner_src  = corner_src;
G.w_min_calib = w_min_calib;
G.cp          = cp;
for j = 1:nP
    G.(PANELS{j, 1}) = Y(:, j);
end

outfile = sprintf('pol1_sweep_%s.mat', SPEC);
save(outfile, '-struct', 'G');
fprintf('Saved %s\n\n', outfile);

end  



function local_vline(x, style, lw)

    yl = ylim();
    ylim(yl);
    plot([x x], yl, style, 'LineWidth', lw);
end

function rep = local_safe_report(x, cfg, cp)

    rep = [];
    try
        r = compute_policy_report_pol(x, cfg, cp);
    catch
        return
    end
    if r.converged
        rep = r;
    end
end

function v = local_pick(rep, fld, scale)
    if isempty(rep)
        v = NaN;
    else
        v = rep.(fld) * scale;
    end
end
