function locate_corner_pol5(SPEC)
% LOCATE_CORNER_POL5  The exact optimum for pol5

if nargin < 1 || isempty(SPEC), SPEC = 'baseline'; end

cfg = pol_config('pol5');
if ~cfg.built
    error('locate_corner_pol5: pol5 is not marked built in pol_config.m.');
end
if ~isequal(cfg.free, {'w_min'})
    error(['locate_corner_pol5: pol_config(''pol5'').free = %s, expected ' ...
           '{w_min} only. '], strjoin(cfg.free, ', '));
end

[cal_file, ~, out_file] = pol_filenames(cfg, SPEC);

% Load the calibration. Nothing here writes to it.
if ~exist(cal_file, 'file')
    error('locate_corner_pol5: %s not found.', cal_file);
end
S = load(cal_file);
if ~isfield(S, 'zeta_star') || ~isfield(S, 'fp')
    error('locate_corner_pol5: %s must contain both ''zeta_star'' and ''fp''.', cal_file);
end
if isfield(S, 'SPEC') && ~isempty(SPEC) && ~strcmp(S.SPEC, SPEC)
    error(['locate_corner_pol5: SPEC MISMATCH. Requested ''%s'' but %s ' ...
           'identifies itself as ''%s''. '], SPEC, cal_file, S.SPEC);
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
cp.tau_R = 0;

fprintf('\n===========================================================\n');
fprintf('  LOCATE_CORNER_POL5  SPEC=''%s''  (tau_R=0, tau_K/tau_L FIXED)\n', SPEC);
fprintf('===========================================================\n');
fprintf('  alpha=%.6f  psi=%.6f  lam_bar=%.6f  eta_bar=%.6f  delta=%.6f\n', ...
        cp.alpha, cp.psi, cp.lam_bar, cp.eta_bar, cp.delta);
fprintf('  r (FIXED) = %.9f   w_min_calib = %.9f\n\n', cp.r, w_min_calib);


welf0 = eval_policy_pol(0, cfg, cp);
if ~welf0.converged || ~strcmp(welf0.regime, 'competitive')
    error(['locate_corner_pol5: evaluating at w_min = 0 did not return a ' ...
           'converged competitive regime (got regime=''%s''. '], ...
          welf0.regime, welf0.converged);
end
threshold  = welf0.w_a1;
SWF_corner = welf0.SWF;

fprintf('--- Threshold wage (competitive w_a1 at w_min = 0) ---\n');
fprintf('  threshold = %.12f   (w_min_calib is %s it)\n', threshold, ...
        tern_str(w_min_calib > threshold, 'ABOVE', 'below'));
fprintf('  regime at w_min=0 = %s   X = %.9f\n', welf0.regime, welf0.X);
fprintf('  (tau_R = 0, so automation is not blocked -- the floor is the only\n');
fprintf('   thing being removed)\n\n');


N_SUBGRID = 41;
w_grid = linspace(0, threshold * (1 - 1e-6), N_SUBGRID);
SWF_w  = -inf(size(w_grid));
reg_w  = repmat({'?'}, size(w_grid));
ws = warning('off', 'all');   

for i = 1:numel(w_grid)
    ww = eval_policy_pol(w_grid(i), cfg, cp);
    reg_w{i} = ww.regime;
    if ww.converged, SWF_w(i) = ww.SWF; end
end
warning(ws);
dev_w = max(abs(SWF_w - SWF_w(1)));
n_not_competitive = sum(~strcmp(reg_w, 'competitive'));

fprintf('--- Flatness check: SWF vs w_min on [0, threshold), %d points ---\n', N_SUBGRID);
fprintf('  max|SWF-SWF(0)| = %.3e   %s\n', dev_w, flat_verdict(dev_w));
if n_not_competitive > 0
    fprintf('  **** %d/%d grid points did NOT return regime=competitive -- unexpected\n', ...
            n_not_competitive, N_SUBGRID);
end
fprintf('\n');


rep_sq = eval_policy_pol(w_min_calib, cfg, cp);
SWF_sq = NaN;
if rep_sq.converged, SWF_sq = rep_sq.SWF; end
fprintf('--- Status quo (floor at its calibrated value) ---\n');
fprintf('  regime = %s   SWF = %.12f   u_a1 = %.9f   w_a1 = %.9f\n\n', ...
        rep_sq.regime, SWF_sq, rep_sq.u_a1, rep_sq.w_a1);


w_min_star = 0;
x_star     = w_min_star;
fprintf('--- Exact optimum ---\n');
fprintf('  w_min* = %.12f (focal point of [0, threshold); the LEVEL is not\n', w_min_star);
fprintf('           identified -- welfare is flat across the whole interval)\n');
fprintf('  SWF*   = %.12f\n', SWF_corner);
fprintf('  u_a1   = %.9f   w_a1 = %.9f   y = %.9f   X = %.9f\n\n', ...
        welf0.u_a1, welf0.w_a1, welf0.y, welf0.X);

save(out_file, 'SPEC', 'w_min_star', 'x_star', 'threshold', 'w_min_calib', ...
     'SWF_corner', 'SWF_sq', 'w_grid', 'SWF_w', 'reg_w', 'dev_w', ...
     'n_not_competitive', 'cp');
fprintf('Saved %s\n\n', out_file);

end  






function s = flat_verdict(dev)
    if dev < 1e-9
        s = 'FLAT to machine precision (plateau confirmed)';
    elseif dev < 1e-6
        s = 'flat to ~1e-6 (plateau confirmed)';
    else
        s = '**** NOT flat -- inspect';
    end
end

function s = tern_str(c, a, b)
    if c, s = a; else, s = b; end
end
