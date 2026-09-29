function locate_corner_pol2(SPEC)
% LOCATE_CORNER_POL2  The regime-flip threshold in closed form

if nargin < 1 || isempty(SPEC), SPEC = 'baseline'; end

cfg = pol_config('pol2');
if ~cfg.built
    error('locate_corner_pol2: pol2 is not marked built in pol_config.m.');
end
if ~isequal(cfg.free, {'tau_K'})
    error(['locate_corner_pol2: pol_config(''pol2'').free = %s, expected ' ...
           '{tau_K} only.'], strjoin(cfg.free, ', '));
end

[cal_file, res_file, out_file] = pol_filenames(cfg, SPEC);



if ~exist(cal_file, 'file')
    error('locate_corner_pol2: %s not found.', cal_file);
end
S = load(cal_file);
if ~isfield(S, 'zeta_star') || ~isfield(S, 'fp')
    error('locate_corner_pol2: %s must contain both ''zeta_star'' and ''fp''.', cal_file);
end
if isfield(S, 'SPEC') && ~isempty(SPEC) && ~strcmp(S.SPEC, SPEC)
    error(['locate_corner_pol2: SPEC MISMATCH. Requested ''%s'' but %s ' ...
           'identifies itself as ''%s''. Refusing to proceed on a ' ...
           'mismatched calibration.'], SPEC, cal_file, S.SPEC);
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
cp.tau_L = fp.tau_L;
cp.w_min = w_min_calib;

alpha   = cp.alpha;
delta   = cp.delta;

[Phi, ~] = phi_from_calibration(fp, false);
cp.Phi = Phi;
tau_K_base = fp.tau_K;

fprintf('\n===========================================================\n');
fprintf('  LOCATE_CORNER_POL2  SPEC=''%s''  (tau_R=0, tau_L, w_min FIXED)\n', SPEC);
fprintf('===========================================================\n');
fprintf('  alpha=%.6f  psi=%.6f  lam_bar=%.6f  eta_bar=%.6f  delta=%.6f\n', ...
        cp.alpha, cp.psi, cp.lam_bar, cp.eta_bar, delta);
fprintf('  w_min (FIXED) = %.9f   tau_L (FIXED) = %.9f   tau_K_base = %.9f   Phi = %.9f\n\n', ...
        w_min_calib, cp.tau_L, tau_K_base, Phi);


r_base_recon = euler_closure(Phi, tau_K_base);
gapA = abs(r_base_recon - fp.r);
fprintf('--- Euler round-trip r(tau_K_base) == fp.r ---\n');
fprintf('  |gap| = %.3e  %s\n\n', gapA, tern(gapA <= 1e-9, 'OK', '**** DRIFT'));



r_base_for_C = euler_closure(Phi, tau_K_base);
[wa1_base, C_star] = noauto_wage_pol(r_base_for_C, cp);
fprintf('--- C (invariant scale constant), probed at r_base ---\n');
fprintf('  w_a1^{c,NA}(base) = %.12f    C = %.12f\n\n', wa1_base, C_star);
fprintf(['  Note: w_a1^{c,NA}(base) %s w_min_calib=%.9f -> regime at ' ...
        'tau_K_base is %s.\n\n'], tern(wa1_base < w_min_calib, '<', '>='), ...
        w_min_calib, tern(wa1_base < w_min_calib, 'mw_binding', 'competitive'));

% =========================================================================
% The flip threshold, in closed form
% =========================================================================
r_flip = (C_star / w_min_calib)^((1 - alpha) / alpha) - delta;
tau_K_flip = 1 - Phi / r_flip;

fprintf('--- tau_K^flip (exact closed form) ---\n');
fprintf('  r_flip = %.9f   tau_K_flip = %.9f\n', r_flip, tau_K_flip);


r_check = euler_closure(Phi, tau_K_flip);
gap_r = abs(r_check - r_flip);
wa1_check = C_star * (r_check + delta)^(-alpha / (1 - alpha));
gap_w = abs(wa1_check - w_min_calib);
fprintf('  round-trip |r_check-r_flip| = %.3e   |w_a1^cNA(flip)-w_min| = %.3e\n', ...
        gap_r, gap_w);
fprintf(['  tau_K < tau_K_flip  ->  regime = competitive (floor slack)\n' ...
         '  tau_K > tau_K_flip  ->  regime = mw_binding  (floor binds)\n' ...
         '  tau_K_flip < 0: a capital SUBSIDY, not a tax (as expected).\n\n']);



TAU_K_LB = -0.50; TAU_K_UB = 1.00;   

N_COARSE = 121;
tk_grid = linspace(TAU_K_LB, TAU_K_UB, N_COARSE);
SWF_grid = -inf(1, N_COARSE);
reg_grid = repmat({'?'}, 1, N_COARSE);
conv_grid = false(1, N_COARSE);
for i = 1:N_COARSE
    try
        ww = eval_policy_pol(tk_grid(i), cfg, cp);
        reg_grid{i} = ww.regime;
        conv_grid(i) = ww.converged;
        if ww.converged, SWF_grid(i) = ww.SWF; end
    catch err
        
        
        reg_grid{i} = 'automation_regime_uncaught_error';
        conv_grid(i) = false;
    end
end
n_conv = sum(conv_grid);
n_conv_below = sum(conv_grid & (tk_grid < tau_K_flip));
n_conv_above = sum(conv_grid & (tk_grid >= tau_K_flip));
n_below = sum(tk_grid < tau_K_flip);
n_above = sum(tk_grid >= tau_K_flip);
fprintf('--- Coarse grid scan (%d points, tau_K in [%.2f, %.2f]) ---\n', ...
        N_COARSE, TAU_K_LB, TAU_K_UB);
fprintf(['  %d/%d points converged in this sandbox (non-convergent points are ' ...
        'automation-regime, which needs fsolve+optimoptions -- MATLAB-only).\n'], ...
        n_conv, N_COARSE);
fprintf(['  Coverage by side of flip: below flip (competitive-eligible) ' ...
        '%d/%d converged; above flip (mw_binding-eligible) %d/%d converged.\n'], ...
        n_conv_below, n_below, n_conv_above, n_above);
if n_conv_below < 0.5 * n_below
    fprintf(['  **** WARNING: coverage below the flip is sparse -- the ' ...
             '"best grid point" below may be UNRELIABLE, not a fair ' ...
             'comparison against the (better-covered) mw_binding side.']);
end
[SWF_best, i_best] = max(SWF_grid);
fprintf('  best grid point: tau_K=%+.6f  SWF=%.10f  regime=%s\n\n', ...
        tk_grid(i_best), SWF_best, reg_grid{i_best});


step = tk_grid(2) - tk_grid(1);
tk_lo = max(TAU_K_LB, tk_grid(i_best) - step);
tk_hi = min(TAU_K_UB, tk_grid(i_best) + step);
tk_ref = linspace(tk_lo, tk_hi, 81);
if tau_K_flip > tk_lo && tau_K_flip < tk_hi
    tk_ref = unique([tk_ref, tau_K_flip - 1e-7, tau_K_flip, tau_K_flip + 1e-7]);
end
SWF_ref = -inf(size(tk_ref)); reg_ref = repmat({'?'}, size(tk_ref));
for i = 1:numel(tk_ref)
    try
        ww = eval_policy_pol(tk_ref(i), cfg, cp);
        reg_ref{i} = ww.regime;
        if ww.converged, SWF_ref(i) = ww.SWF; end
    catch
        reg_ref{i} = 'automation_regime_uncaught_error';
    end
end
[SWF_star_grid, j_best] = max(SWF_ref);
tau_K_star_grid = tk_ref(j_best);
fprintf('--- Refined grid candidate ---\n');
fprintf('  tau_K*(grid) = %+.9f   SWF*(grid) = %.12f   regime=%s\n', ...
        tau_K_star_grid, SWF_star_grid, reg_ref{j_best});
if abs(tau_K_star_grid - tau_K_flip) < 2*step
    fprintf(['  >> grid candidate sits NEAR the flip -- consistent with, but ' ...
            'NOT PROOF of, the optimum sitting at or just past the flip. ' ...
            'This is a %d-point grid, not a refined search.\n'], numel(tk_ref));
end
fprintf(['  **** THIS IS A GRID-QUANTIZED CANDIDATE, NOT A VERIFIED OPTIMUM. ****\n' ...
         '  Run welfare_search_pol(''pol2'', SPEC)''s Phase 2 SA (MATLAB-only) ' ...
         'for the real answer to whether tau_K* sits strictly beyond ' ...
         'tau_K_flip or is bounded by grid/search resolution near it -- ' ...
         'this is an open empirical question.\n\n']);


sa_tau_K = NaN; sa_swf = NaN;
if exist(res_file, 'file')
    R2 = load(res_file);
    if all(isfield(R2, {'tau_K_star','SWF_star'}))
        sa_tau_K = R2.tau_K_star; sa_swf = R2.SWF_star;
        fprintf('--- SA cross-check ---\n');
        fprintf('  SA: tau_K*=%.6f  SWF*=%.10f\n', sa_tau_K, sa_swf);
        fprintf('  gap vs. grid candidate: tau_K %.3e, SWF %.3e\n\n', ...
                abs(sa_tau_K - tau_K_star_grid), abs(sa_swf - SWF_star_grid));
    end
end

is_feasibility_boundary_only = true;
is_grid_quantized_only = true;
save(out_file, 'SPEC', 'Phi', 'C_star', 'tau_K_flip', 'r_flip', 'gap_r', 'gap_w', ...
     'tk_grid', 'SWF_grid', 'reg_grid', 'conv_grid', ...
     'n_conv_below', 'n_below', 'n_conv_above', 'n_above', ...
     'tk_ref', 'SWF_ref', 'reg_ref', 'tau_K_star_grid', 'SWF_star_grid', ...
     'w_min_calib', 'sa_tau_K', 'sa_swf', 'cp', ...
     'is_feasibility_boundary_only', 'is_grid_quantized_only');
fprintf('Saved %s\n\n', out_file);

end  




function s = tern(cond, s_true, s_false)
    if cond, s = s_true; else, s = s_false; end
end
