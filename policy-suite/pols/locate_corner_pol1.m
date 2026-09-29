function locate_corner_pol1(SPEC)
% LOCATE_CORNER_POL1  The automation-blocking corner for pol1

%
%   locate_corner_pol1('baseline')

if nargin < 1 || isempty(SPEC), SPEC = 'baseline'; end

cfg = pol_config('pol1');
if ~cfg.built
    error('locate_corner_pol1: pol1 is not marked built in pol_config.m.');
end
if ~isequal(cfg.free, {'tau_R'})
    error(['locate_corner_pol1: pol_config(''pol1'').free = %s, expected ' ...
           '{tau_R} only. '], ...
          strjoin(cfg.free, ', '));
end

[cal_file, res_file, out_file] = pol_filenames(cfg, SPEC);


if ~exist(cal_file, 'file')
    error('locate_corner_pol1: %s not found.', cal_file);
end
S = load(cal_file);
if ~isfield(S, 'zeta_star') || ~isfield(S, 'fp')
    error('locate_corner_pol1: %s must contain both ''zeta_star'' and ''fp''.', cal_file);
end
if isfield(S, 'SPEC') && ~isempty(SPEC) && ~strcmp(S.SPEC, SPEC)
    error(['locate_corner_pol1: SPEC MISMATCH. Requested ''%s'' but %s ' ...
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

alpha   = cp.alpha;
delta   = cp.delta;
lam_bar = cp.lam_bar;
rd      = cp.r + delta;

fprintf('\n===========================================================\n');
fprintf('  LOCATE_CORNER_POL1  SPEC=''%s''  (tau_K, tau_L, w_min ALL FIXED)\n', SPEC);
fprintf('===========================================================\n');
fprintf('  alpha=%.6f  psi=%.6f  lam_bar=%.6f  eta_bar=%.6f  delta=%.6f\n', ...
        cp.alpha, cp.psi, cp.lam_bar, cp.eta_bar, delta);
fprintf('  r (FIXED) = %.9f   r+delta = %.9f   w_min (FIXED) = %.9f\n\n', ...
        cp.r, rd, w_min_calib);


[wa1_cNA, C_star] = noauto_wage_pol(cp.r, cp);
fprintf('--- w_a1^{c,NA} (no-automation competitive wage at r) ---\n');
fprintf('  w_a1^{c,NA} = %.12f    C = %.12f\n\n', wa1_cNA, C_star);

[corner, branch, wa1_eff] = corner_of_pol(cp.r, delta, lam_bar, w_min_calib, C_star, alpha);
fprintf('--- tau_R^corner at w_min = w_min_calib ---\n');
fprintf('  branch=%s   w_a1_eff=%.9f   tau_R^corner = %.9f\n\n', ...
        branch, wa1_eff, corner);
if strcmp(branch, 'binding')
    fprintf(['  NOTE: branch=binding (w_min_calib > w_a1^{c,NA}): pol1 asks "block automation ' ...
             'with tau_R alone, GIVEN the wage floor is where it actually ' ...
             'is," while pol6 asks "block automation with tau_R AND the ' ...
             'wage floor jointly chosen."\n\n']);
end



PLATEAU_EPS = 1e-6;
tau_grid = corner + [PLATEAU_EPS, 1e-3, 0.05, 0.25, 1.0, 2.5];
SWF_t = -inf(size(tau_grid));
for i = 1:numel(tau_grid)
    x = tau_grid(i);   % single free instrument
    ww = eval_policy_pol(x, cfg, cp);
    if ww.converged, SWF_t(i) = ww.SWF; end
end
dev_t = max(abs(SWF_t - SWF_t(1)));
fprintf('--- Flatness check: SWF vs tau_R above corner (w_min fixed) ---\n');
fprintf('  max|SWF-SWF(corner)| = %.3e   %s\n\n', dev_t, flat_verdict(dev_t));

SWF_corner = SWF_t(1);
fprintf('--- Exact optimum ---\n');
fprintf('  tau_R* = %.12f   w_min (FIXED) = %.9f   SWF* = %.12f\n\n', ...
        corner, w_min_calib, SWF_corner);


sa_tau = NaN; sa_swf = NaN;
if exist(res_file, 'file')
    R2 = load(res_file);
    if all(isfield(R2, {'tau_R_star','SWF_star'}))
        sa_tau = R2.tau_R_star; sa_swf = R2.SWF_star;
        fprintf('--- SA cross-check ---\n');
        fprintf('  SA: tau_R*=%.6f  SWF*=%.10f\n', sa_tau, sa_swf);
        fprintf('  |SWF diff| = %.3e\n\n', abs(sa_swf - SWF_corner));
    end
end

save(out_file, 'SPEC', 'corner', 'branch', 'wa1_eff', 'wa1_cNA', 'C_star', ...
     'w_min_calib', 'SWF_corner', 'tau_grid', 'SWF_t', 'dev_t', ...
     'sa_tau', 'sa_swf', 'cp');
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
