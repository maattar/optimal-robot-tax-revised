function locate_corner_pol(pol_id, SPEC, tau_L_override)
% LOCATE_CORNER_POL  The exact automation-blocking corner for pol6 and
% pol4, taking its shape from pol_config.



if nargin < 2 || isempty(SPEC),           SPEC = 'baseline'; end
if nargin < 3,                            tau_L_override = []; end

cfg = pol_config(pol_id);
if ~cfg.built
    error(['locate_corner_pol: %s is not marked built in pol_config.'], pol_id);
end
if ~any(strcmp(pol_id, {'pol6', 'pol4'}))
    error(['locate_corner_pol: %s is not supported here.'], pol_id);
end

[cal_file, res_file, out_file] = pol_filenames(cfg, SPEC);


if ~exist(cal_file, 'file')
    error('locate_corner_pol: %s not found.', cal_file);
end
S = load(cal_file);
if ~isfield(S, 'zeta_star') || ~isfield(S, 'fp')
    error('locate_corner_pol: %s must contain both ''zeta_star'' and ''fp''.', cal_file);
end
if isfield(S, 'SPEC') && ~isempty(SPEC) && ~strcmp(S.SPEC, SPEC)
    error(['locate_corner_pol: SPEC MISMATCH. Requested ''%s'' but %s ' ...
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

alpha   = cp.alpha;
delta   = cp.delta;
lam_bar = cp.lam_bar;

fprintf('\n===========================================================\n');
fprintf('  LOCATE_CORNER_POL(''%s'')  SPEC=''%s''\n', pol_id, SPEC);
fprintf('===========================================================\n');
fprintf('  alpha=%.6f  psi=%.6f  lam_bar=%.6f  eta_bar=%.6f  delta=%.6f\n', ...
        cp.alpha, cp.psi, cp.lam_bar, cp.eta_bar, delta);

if cfg.needs_euler
    
    
    [Phi, ~] = phi_from_calibration(fp, false);
    cp.Phi   = Phi;
    cp.w_min = w_min_calib;
    tau_K_base = fp.tau_K;

    fprintf('  w_min (FIXED) = %.9f   tau_K_base = %.9f   Phi = %.9f\n\n', ...
            w_min_calib, tau_K_base, Phi);

    
    r_base_recon = euler_closure(Phi, tau_K_base);
    gapA = abs(r_base_recon - fp.r);
    fprintf('--- Euler round-trip r(tau_K_base) == fp.r ---\n');
    fprintf('  |gap| = %.3e  %s\n\n', gapA, tern(gapA <= 1e-9, 'OK', '**** DRIFT'));

    
    r_base_for_C = euler_closure(Phi, tau_K_base);
    [wa1_base, C_star] = noauto_wage_pol(r_base_for_C, cp);
    fprintf('--- C (invariant scale constant), probed at r_base ---\n');
    fprintf('  w_a1^{c,NA}(base) = %.12f    C = %.12f\n\n', wa1_base, C_star);

   
    
    if ~isempty(tau_L_override)
        tau_L_use = tau_L_override; tau_L_src = 'explicit override';
    elseif exist(res_file, 'file')
        R = load(res_file);
        if isfield(R, 'tau_L_star')
            tau_L_use = R.tau_L_star; tau_L_src = sprintf('tau_L* from %s', res_file);
        else
            tau_L_use = fp.tau_L; tau_L_src = 'fp.tau_L (result file lacks tau_L_star)';
        end
    else
        tau_L_use = fp.tau_L; tau_L_src = 'fp.tau_L (no search result found)';
        warning(['locate_corner_pol: no welfare_search result found']);
    end
    fprintf('  tau_L held at %.9f  (%s)\n\n', tau_L_use, tau_L_src);


    r_base = euler_closure(Phi, tau_K_base);
    [corner_base, branch_base] = corner_of_pol(r_base, delta, lam_bar, w_min_calib, C_star, alpha);
    fprintf('--- Corner at tau_K_base (should equal the tau_R-only suite''s corner) ---\n');
    fprintf('  branch=%s   tau_R^corner(base) = %.9f\n\n', branch_base, corner_base);


    
    rd_flip = (C_star / w_min_calib)^((1 - alpha) / alpha);
    r_flip  = rd_flip - delta;
    if r_flip > 0
        tau_K_flip = 1 - Phi / r_flip;
        fprintf('  tau_K^flip (branch switch) = %+.9f\n\n', tau_K_flip);
    else
        tau_K_flip = NaN;
        fprintf('  No flip: floor binds for every admissible tau_K.\n\n');
    end

    
    TAU_K_LB = -0.5; TAU_K_UB = 1.0; TAU_R_LB = -0.10; PLATEAU_EPS = 1e-6;
    N_CURVE = 61;
    tk_curve  = linspace(TAU_K_LB, TAU_K_UB, N_CURVE);
    SWF_curve = -inf(1, N_CURVE);
    corner_curve = nan(1, N_CURVE);
    reg_curve = repmat({'?'}, 1, N_CURVE);
    for i = 1:N_CURVE
        r_i = euler_closure(Phi, tk_curve(i));
        corner_curve(i) = corner_of_pol(r_i, delta, lam_bar, w_min_calib, C_star, alpha);
        x = [max(corner_curve(i), TAU_R_LB) + PLATEAU_EPS, tk_curve(i), tau_L_use];
        ww = eval_policy_pol(x, cfg, cp);
        reg_curve{i} = ww.regime;
        if ww.converged, SWF_curve(i) = ww.SWF; end
    end
    [SWF_curve_max, i_max] = max(SWF_curve);
    fprintf('--- Frontier scan (%d points, tau_K in [%.2f, %.2f]) ---\n', ...
            N_CURVE, TAU_K_LB, TAU_K_UB);
    fprintf('  best grid point: tau_K=%+.6f  tau_R^corner=%.6f  SWF=%.10f  regime=%s\n\n', ...
            tk_curve(i_max), corner_curve(i_max), SWF_curve_max, reg_curve{i_max});
 
    step_K = tk_curve(2) - tk_curve(1);
    tk_lo  = max(TAU_K_LB, tk_curve(i_max) - step_K);
    tk_hi  = min(TAU_K_UB, tk_curve(i_max) + step_K);
    tk_ref = linspace(tk_lo, tk_hi, 81);
    if ~isnan(tau_K_flip) && tau_K_flip > tk_lo && tau_K_flip < tk_hi
        tk_ref = unique([tk_ref, tau_K_flip - 1e-7, tau_K_flip, tau_K_flip + 1e-7]);
    end
    SWF_ref = -inf(size(tk_ref)); cr_ref = nan(size(tk_ref));
    reg_ref = repmat({'?'}, size(tk_ref));
    for i = 1:numel(tk_ref)
        r_i = euler_closure(Phi, tk_ref(i));
        cr_ref(i) = corner_of_pol(r_i, delta, lam_bar, w_min_calib, C_star, alpha);
        x = [max(cr_ref(i), TAU_R_LB) + PLATEAU_EPS, tk_ref(i), tau_L_use];
        ww = eval_policy_pol(x, cfg, cp);
        reg_ref{i} = ww.regime;
        if ww.converged, SWF_ref(i) = ww.SWF; end
    end
    [SWF_star, j_max] = max(SWF_ref);
    tau_K_star = tk_ref(j_max);
    tau_R_star = cr_ref(j_max);
    fprintf('--- Refined optimum ---\n');
    fprintf('  tau_K* = %+.9f   tau_R* = %.9f   SWF* = %.12f   regime=%s\n', ...
            tau_K_star, tau_R_star, SWF_star, reg_ref{j_max});
    if ~isnan(tau_K_flip) && abs(tau_K_star - tau_K_flip) < 1e-6
        fprintf('  >> optimum sits ON the regime flip (kink, not a smooth interior max)\n');
    end
    fprintf('\n');

    
    tau_grid = corner_base + [0, 1e-3, 0.05, 0.25, 1.0, 2.5];
    SWF_t = -inf(size(tau_grid));
    for i = 1:numel(tau_grid)
        x = [tau_grid(i), tau_K_base, tau_L_use];
        ww = eval_policy_pol(x, cfg, cp);
        if ww.converged, SWF_t(i) = ww.SWF; end
    end
    dev_t = max(abs(SWF_t - SWF_t(1)));
    fprintf('--- Flatness check: SWF vs tau_R above corner (tau_K_base) ---\n');
    fprintf('  max|SWF-SWF(corner)| = %.3e   %s\n\n', dev_t, flat_verdict(dev_t));

    
    sa_tau_R = NaN; sa_tau_K = NaN; sa_tau_L = NaN; sa_swf = NaN;
    if exist(res_file, 'file')
        R2 = load(res_file);
        if all(isfield(R2, {'tau_R_star','tau_K_star','tau_L_star','SWF_star'}))
            sa_tau_R = R2.tau_R_star; sa_tau_K = R2.tau_K_star;
            sa_tau_L = R2.tau_L_star; sa_swf = R2.SWF_star;
            fprintf('--- SA cross-check ---\n');
            fprintf('  SA: tau_R=%.6f tau_K=%+.6f tau_L=%+.6f SWF=%.10f\n', ...
                    sa_tau_R, sa_tau_K, sa_tau_L, sa_swf);
            fprintf('  gap tau_R*(SA)-corner(SA''s tau_K) reported in the .mat only.\n\n');
        end
    end

    save(out_file, 'SPEC', 'Phi', 'C_star', 'corner_base', 'branch_base', ...
         'tau_K_flip', 'tk_curve', 'corner_curve', 'SWF_curve', 'reg_curve', ...
         'tk_ref', 'cr_ref', 'SWF_ref', 'reg_ref', ...
         'tau_K_star', 'tau_R_star', 'SWF_star', 'tau_L_use', 'tau_L_src', ...
         'w_min_calib', 'sa_tau_R', 'sa_tau_K', 'sa_tau_L', 'sa_swf', 'cp');
    fprintf('Saved %s\n\n', out_file);

else
    
    
    cp.r     = fp.r;
    cp.tau_K = fp.tau_K;
    cp.tau_L = fp.tau_L;
    rd = cp.r + delta;

    fprintf('  r (FIXED) = %.9f   r+delta = %.9f\n\n', cp.r, rd);

    [wa1_cNA, C_star] = noauto_wage_pol(cp.r, cp);
    fprintf('--- w_a1^{c,NA} and corner at w_min = 0 ---\n');
    fprintf('  w_a1^{c,NA} = %.12f\n', wa1_cNA);

    [corner, branch] = corner_of_pol(cp.r, delta, lam_bar, 0, C_star, alpha);
    fprintf('  branch=%s   tau_R^corner = %.9f\n\n', branch, corner);

    
    
    tR_plateau = corner + max(0.05, 0.05 * abs(corner));
    w_grid = linspace(0, wa1_cNA, 11);
    SWF_w = -inf(size(w_grid));
    for i = 1:numel(w_grid)
        x = [tR_plateau, w_grid(i)];
        ww = eval_policy_pol(x, cfg, cp);
        if ww.converged, SWF_w(i) = ww.SWF; end
    end
    dev_w = max(abs(SWF_w - SWF_w(1)));
    fprintf('--- Flatness check: SWF vs w_min on [0, w_a1^c] ---\n');
    fprintf('  max|SWF-SWF(w_min=0)| = %.3e   %s\n\n', dev_w, flat_verdict(dev_w));

    
    
    PLATEAU_EPS = 1e-6;
    tau_grid = corner + [PLATEAU_EPS, 1e-3, 0.05, 0.25, 1.0, 2.5];
    SWF_t = -inf(size(tau_grid));
    for i = 1:numel(tau_grid)
        x = [tau_grid(i), 0];
        ww = eval_policy_pol(x, cfg, cp);
        if ww.converged, SWF_t(i) = ww.SWF; end
    end
    dev_t = max(abs(SWF_t - SWF_t(1)));
    fprintf('--- Flatness check: SWF vs tau_R above corner (w_min=0) ---\n');
    fprintf('  max|SWF-SWF(corner)| = %.3e   %s\n\n', dev_t, flat_verdict(dev_t));

    SWF_corner = SWF_t(1);
    fprintf('--- Exact optimum ---\n');
    fprintf('  tau_R* = %.12f   w_min* = 0   SWF* = %.12f\n\n', corner, SWF_corner);

    sa_tau = NaN; sa_w = NaN; sa_swf = NaN;
    if exist(res_file, 'file')
        R2 = load(res_file);
        if all(isfield(R2, {'tau_R_star','w_min_star','SWF_star'}))
            sa_tau = R2.tau_R_star; sa_w = R2.w_min_star; sa_swf = R2.SWF_star;
            fprintf('--- SA cross-check ---\n');
            fprintf('  SA: tau_R*=%.6f  w_min*=%.6f  SWF*=%.10f\n', sa_tau, sa_w, sa_swf);
            fprintf('  |SWF diff| = %.3e\n\n', abs(sa_swf - SWF_corner));
        end
    end

    save(out_file, 'SPEC', 'corner', 'branch', 'wa1_cNA', 'C_star', ...
         'SWF_corner', 'w_grid', 'SWF_w', 'dev_w', 'tau_grid', 'SWF_t', 'dev_t', ...
         'sa_tau', 'sa_w', 'sa_swf', 'cp');
    fprintf('Saved %s\n\n', out_file);
end

end  



function s = tern(cond, s_true, s_false)
    if cond, s = s_true; else, s = s_false; end
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
