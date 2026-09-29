function locate_corner_pol3(SPEC)
% LOCATE_CORNER_POL3  pol2's regime-flip threshold


if nargin < 1 || isempty(SPEC), SPEC = 'baseline'; end

cfg = pol_config('pol3');
if ~cfg.built
    error('locate_corner_pol3: pol3 is not marked built in pol_config.m.');
end
if ~isequal(cfg.free, {'tau_K', 'tau_L'})
    error(['locate_corner_pol3: pol_config(''pol3'').free = %s, expected ' ...
           '{tau_K, tau_L}.'], strjoin(cfg.free, ', '));
end

[cal_file, res_file, out_file] = pol_filenames(cfg, SPEC);


if ~exist(cal_file, 'file')
    error('locate_corner_pol3: %s not found.', cal_file);
end
S = load(cal_file);
if ~isfield(S, 'zeta_star') || ~isfield(S, 'fp')
    error('locate_corner_pol3: %s must contain both ''zeta_star'' and ''fp''.', cal_file);
end
if isfield(S, 'SPEC') && ~isempty(SPEC) && ~strcmp(S.SPEC, SPEC)
    error(['locate_corner_pol3: SPEC MISMATCH. Requested ''%s'' but %s ' ...
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


cp.tau_R = 0;
cp.w_min = w_min_calib;

alpha   = cp.alpha;
delta   = cp.delta;

[Phi, ~] = phi_from_calibration(fp, false);
cp.Phi = Phi;
tau_K_base = fp.tau_K;

fprintf('\n===========================================================\n');
fprintf('  LOCATE_CORNER_POL3  SPEC=''%s''  (tau_R=0, w_min=%.9f FIXED; tau_K, tau_L FREE)\n', SPEC, w_min_calib);
fprintf('===========================================================\n');
fprintf('  alpha=%.6f  psi=%.6f  lam_bar=%.6f  eta_bar=%.6f  delta=%.6f\n', ...
        cp.alpha, cp.psi, cp.lam_bar, cp.eta_bar, delta);
fprintf('  w_min (FIXED, status quo) = %.9f   tau_L_base = %.9f   tau_K_base = %.9f   Phi = %.9f\n\n', ...
        w_min_calib, fp.tau_L, tau_K_base, Phi);


r_base_recon = euler_closure(Phi, tau_K_base);
gapA = abs(r_base_recon - fp.r);
fprintf('--- Euler round-trip r(tau_K_base) == fp.r ---\n');
fprintf('  |gap| = %.3e  %s\n\n', gapA, tern(gapA <= 1e-9, 'OK', '**** DRIFT'));


[wa1_base, C_star] = noauto_wage_pol(r_base_recon, cp);
fprintf('--- C (invariant scale constant), probed at r_base ---\n');
fprintf('  w_a1^{c,NA}(base) = %.12f    C = %.12f\n\n', wa1_base, C_star);
fprintf(['  Note: w_a1^{c,NA}(base) %s w_min_calib=%.9f -> regime at ' ...
        'tau_K_base is %s (tau_L never enters this comparison).\n\n'], ...
        tern(wa1_base < w_min_calib, '<', '>='), w_min_calib, ...
        tern(wa1_base < w_min_calib, 'mw_binding', 'competitive'));

% =========================================================================
% THE FLIP THRESHOLD (exact, closed form -- identical to pol2's)
% =========================================================================
r_flip = (C_star / w_min_calib)^((1 - alpha) / alpha) - delta;
tau_K_flip = 1 - Phi / r_flip;

fprintf('--- tau_K^flip (exact closed form, IDENTICAL to pol2''s) ---\n');
fprintf('  r_flip = %.9f   tau_K_flip = %.9f\n', r_flip, tau_K_flip);

r_check = euler_closure(Phi, tau_K_flip);
gap_r = abs(r_check - r_flip);
wa1_check = C_star * (r_check + delta)^(-alpha / (1 - alpha));
gap_w = abs(wa1_check - w_min_calib);
fprintf('  round-trip |r_check-r_flip| = %.3e   |w_a1^cNA(flip)-w_min| = %.3e\n', ...
        gap_r, gap_w);
fprintf(['  tau_K < tau_K_flip  ->  regime = competitive (floor slack)\n' ...
         '  tau_K > tau_K_flip  ->  regime = mw_binding  (floor binds)\n\n']);


[~, ~, out_file2] = pol_filenames(pol_config('pol2'), SPEC);
if exist(out_file2, 'file')
    C2 = load(out_file2);
    if isfield(C2, 'tau_K_flip')
        gapK = abs(C2.tau_K_flip - tau_K_flip);
        fprintf('--- Cross-suite check vs. pol2''s tau_K_flip (%s) ---\n', out_file2);
        fprintf('  pol2 tau_K_flip = %.9f   pol3 tau_K_flip = %.9f   |gap| = %.3e  %s\n\n', ...
                C2.tau_K_flip, tau_K_flip, gapK, ...
                tern(gapK <= 1e-9, 'MATCH (as required)', '**** MISMATCH -- INVESTIGATE'));
    end
else
    fprintf(['--- Cross-suite check vs. pol2 skipped: %s not found (run ' ...
             'locate_corner_pol2(''%s'') to enable it) ---\n\n'], out_file2, SPEC);
end




PLATEAU_EPS = 1e-6;
tauL_probe = [fp.tau_L - 0.10, fp.tau_L, fp.tau_L + 0.10, 0.0];
fprintf('--- tau_L-independence check of the regime label at tau_K_flip +/- eps ---\n');
regimes_below = cell(size(tauL_probe));
regimes_above = cell(size(tauL_probe));
for i = 1:numel(tauL_probe)
   
    
    try
        w_below = eval_policy_pol([tau_K_flip - PLATEAU_EPS, tauL_probe(i)], cfg, cp);
        regimes_below{i} = w_below.regime;
    catch
        regimes_below{i} = 'automation_regime_uncaught_error';
    end
    try
        w_above = eval_policy_pol([tau_K_flip + PLATEAU_EPS, tauL_probe(i)], cfg, cp);
        regimes_above{i} = w_above.regime;
    catch
        regimes_above{i} = 'automation_regime_uncaught_error';
    end
    fprintf('  tau_L=%+.6f :  regime(tau_K_flip-eps)=%-14s  regime(tau_K_flip+eps)=%s\n', ...
            tauL_probe(i), regimes_below{i}, regimes_above{i});
end
below_ok = all(strcmp(regimes_below, regimes_below{1}));
above_ok = all(strcmp(regimes_above, regimes_above{1}));
fprintf('  %s\n\n', tern(below_ok && above_ok, ...
    'CONFIRMED: regime label invariant to tau_L on both sides of the flip'));


TAU_K_LB = -0.50; TAU_K_UB = 1.00;   % matches welfare_search_pol.m's DEFAULT_BOUNDS.tau_K
TAU_L_LB = -0.50; TAU_L_UB = 1.00;   % matches welfare_search_pol.m's DEFAULT_BOUNDS.tau_L
N_K = 121; N_L = 41;
tk_grid = linspace(TAU_K_LB, TAU_K_UB, N_K);
tl_grid = linspace(TAU_L_LB, TAU_L_UB, N_L);

SWF_grid = -inf(N_K, N_L);
reg_grid = repmat({'?'}, N_K, N_L);
conv_grid = false(N_K, N_L);
for i = 1:N_K
    for j = 1:N_L
        try
            ww = eval_policy_pol([tk_grid(i), tl_grid(j)], cfg, cp);
            reg_grid{i, j} = ww.regime;
            conv_grid(i, j) = ww.converged;
            if ww.converged, SWF_grid(i, j) = ww.SWF; end
        catch
            
            
            reg_grid{i, j} = 'automation_regime_uncaught_error';
            conv_grid(i, j) = false;
        end
    end
end
n_conv = sum(conv_grid(:));
N_TOTAL = N_K * N_L;
fprintf('--- Coarse 2-D grid scan (%d x %d = %d points, tau_K in [%.2f, %.2f], tau_L in [%.2f, %.2f]) ---\n', ...
        N_K, N_L, N_TOTAL, TAU_K_LB, TAU_K_UB, TAU_L_LB, TAU_L_UB);
fprintf(['  %d/%d points converged in this sandbox (non-convergent points are ' ...
        'automation-regime, which needs fsolve+optimoptions -- MATLAB-only).\n\n'], ...
        n_conv, N_TOTAL);

[SWF_best, lin_idx] = max(SWF_grid(:));
[i_best, j_best] = ind2sub(size(SWF_grid), lin_idx);
fprintf('  best grid point: tau_K=%+.6f  tau_L=%+.6f  SWF=%.10f  regime=%s\n\n', ...
        tk_grid(i_best), tl_grid(j_best), SWF_best, reg_grid{i_best, j_best});


step_k = tk_grid(2) - tk_grid(1);
step_l = tl_grid(2) - tl_grid(1);
tk_lo = max(TAU_K_LB, tk_grid(i_best) - step_k);
tk_hi = min(TAU_K_UB, tk_grid(i_best) + step_k);
tl_lo = max(TAU_L_LB, tl_grid(j_best) - step_l);
tl_hi = min(TAU_L_UB, tl_grid(j_best) + step_l);
tk_ref = linspace(tk_lo, tk_hi, 21);
tl_ref = linspace(tl_lo, tl_hi, 21);
if tau_K_flip > tk_lo && tau_K_flip < tk_hi
    tk_ref = unique([tk_ref, tau_K_flip - PLATEAU_EPS, tau_K_flip, tau_K_flip + PLATEAU_EPS]);
end

SWF_ref = -inf(numel(tk_ref), numel(tl_ref));
reg_ref = repmat({'?'}, numel(tk_ref), numel(tl_ref));
for i = 1:numel(tk_ref)
    for j = 1:numel(tl_ref)
        try
            ww = eval_policy_pol([tk_ref(i), tl_ref(j)], cfg, cp);
            reg_ref{i, j} = ww.regime;
            if ww.converged, SWF_ref(i, j) = ww.SWF; end
        catch
            reg_ref{i, j} = 'automation_regime_uncaught_error';
        end
    end
end
[SWF_star_grid, lin_idx2] = max(SWF_ref(:));
[i2, j2] = ind2sub(size(SWF_ref), lin_idx2);
tau_K_star_grid = tk_ref(i2);
tau_L_star_grid = tl_ref(j2);
fprintf('--- Refined grid candidate ---\n');
fprintf('  tau_K*(grid) = %+.9f   tau_L*(grid) = %+.9f   SWF*(grid) = %.12f   regime=%s\n', ...
        tau_K_star_grid, tau_L_star_grid, SWF_star_grid, reg_ref{i2, j2});
fprintf(['  **** THIS IS A GRID-QUANTIZED CANDIDATE, NOT A VERIFIED OPTIMUM. ****\n' ...
         '  Run welfare_search_pol(''pol3'', SPEC)''s Phase 2 SA (MATLAB-only) ' ...
         'for the real answer.\n\n']);


sa_tau_K = NaN; sa_tau_L = NaN; sa_swf = NaN;
if exist(res_file, 'file')
    R2 = load(res_file);
    if all(isfield(R2, {'tau_K_star', 'tau_L_star', 'SWF_star'}))
        sa_tau_K = R2.tau_K_star; sa_tau_L = R2.tau_L_star; sa_swf = R2.SWF_star;
        fprintf('--- SA cross-check ---\n');
        fprintf('  SA: tau_K*=%.6f  tau_L*=%.6f  SWF*=%.10f\n', sa_tau_K, sa_tau_L, sa_swf);
        fprintf('  gap vs. grid candidate: tau_K %.3e, tau_L %.3e, SWF %.3e\n\n', ...
                abs(sa_tau_K - tau_K_star_grid), abs(sa_tau_L - tau_L_star_grid), abs(sa_swf - SWF_star_grid));
    end
end

is_feasibility_boundary_only = true;
is_grid_quantized_only = true;
save(out_file, 'SPEC', 'Phi', 'C_star', 'tau_K_flip', 'r_flip', 'gap_r', 'gap_w', ...
     'tk_grid', 'tl_grid', 'SWF_grid', 'reg_grid', 'conv_grid', ...
     'tk_ref', 'tl_ref', 'SWF_ref', 'reg_ref', ...
     'tau_K_star_grid', 'tau_L_star_grid', 'SWF_star_grid', ...
     'w_min_calib', 'sa_tau_K', 'sa_tau_L', 'sa_swf', 'cp', ...
     'is_feasibility_boundary_only', 'is_grid_quantized_only');
fprintf('Saved %s\n\n', out_file);

end 





function s = tern(cond, s_true, s_false)
    if cond, s = s_true; else, s = s_false; end
end
