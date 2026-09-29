function sol = solve_minwage_bgp(params, w_min)
% SOLVE_MINWAGE_BGP   binds fortype-1 labour.
%
%   sol = solve_minwage_bgp(params, w_min)

sc    = compute_scalars(params);
psi_A = sc.psi * sc.A;       
gam   = sc.gamma;             % theta/(1-theta), between -1 and 0

% automatio
hat_lam_fixed = sc.p / w_min;
auto_on       = (hat_lam_fixed < sc.lam_bar);

if params.verbose
    fprintf('\n=== Min-Wage BGP Solver (w_min = %.6f) ===\n', w_min);
    if auto_on
        fprintf('Exogenous automation check: AUTOMATION\n');
        fprintf('  hat_lambda = %.6f  <  lam_bar = %.6f\n', ...
                hat_lam_fixed, sc.lam_bar);
    else
        fprintf('Exogenous automation check: NO AUTOMATION\n');
        fprintf('  hat_lambda = %.6f  >= lam_bar = %.6f\n', ...
                hat_lam_fixed, sc.lam_bar);
    end
end

% regimes
if auto_on
    regime_list = {'case1', 'case2'};    % the first case, then the second
else
    regime_list = {'no_auto'};            % nothing is automated
end

regime_labels.no_auto = 'no-automation';
regime_labels.case1   = 'Case 1 (ray exits top)';
regime_labels.case2   = 'Case 2 (ray exits right)';

% fzero settings
fzero_opts      = optimset('TolX', 1e-12, 'Display', 'off');
RESID_TOL       = 1e-8;        % the residual must clear this as well
MAX_BRACK_ITERS = 40;          % how far the search for an interval may go


v_lo_base = log(w_min * 1.001);

% regimes
sol_found  = false;
sol_regime = '';
v_sol      = NaN;

for ri = 1:numel(regime_list)
    regime = regime_list{ri};

    if params.verbose
        fprintf('\nTrying regime: %s ...\n', regime);
    end

    res_fun = @(v) residual_minwage(v, w_min, regime, sc);


    bracket_ok = false;
    v_prev = v_lo_base;
    R_prev = res_fun(v_prev);

    for iter = 1:MAX_BRACK_ITERS
        v_next = v_prev + log(2);      % double the wage
        R_next = res_fun(v_next);

        if abs(R_prev) < 1e5 && abs(R_next) < 1e5 ...
                && sign(R_prev) ~= sign(R_next)
            v_lo       = v_prev;
            v_hi       = v_next;
            bracket_ok = true;
            break
        end

        % Move on
        v_prev = v_next;
        R_prev = R_next;
    end

    if ~bracket_ok
        if params.verbose
            fprintf('  Could not bracket residual for regime %s. Skipping.\n', regime);
        end
        continue
    end

    if params.verbose
        fprintf('  Bracket found: log(w_a2) in [%.4f, %.4f]  (w_a2 in [%.4f, %.4f])\n', ...
                v_lo, v_hi, exp(v_lo), exp(v_hi));
    end

    %fzero   call
    try
        [v_sol_try, ~, exitflag] = fzero(res_fun, [v_lo, v_hi], fzero_opts);
    catch ME
        if params.verbose
            fprintf('  fzero error for regime %s: %s\n', regime, ME.message);
        end
        continue
    end

    if exitflag <= 0
        if params.verbose
            fprintf('  fzero did not converge for regime %s (exitflag=%d). Skipping.\n', ...
                    regime, exitflag);
        end
        continue
    end

       %  check
    R_check = res_fun(v_sol_try);
    if abs(R_check) > RESID_TOL
        if params.verbose
            fprintf('  |R| = %.2e > %.2e at fzero solution. Skipping regime %s.\n', ...
                    abs(R_check), RESID_TOL, regime);
        end
        continue
    end

    % consistency-
    w_a2_try = exp(v_sol_try);

    regime_ok = false;
    switch regime
        case 'no_auto'
            
            hat_eta_try = w_a2_try / w_min;
            regime_ok   = (hat_eta_try < sc.eta_bar);
            if params.verbose && ~regime_ok
                fprintf('  Regime no_auto fails: hat_eta = %.4f >= bar_eta = %.4f\n', ...
                        hat_eta_try, sc.eta_bar);
            end

        case 'case1'
              cond1     = (hat_lam_fixed < sc.lam_bar);
            cond2     = (sc.lam_bar < (sc.p / w_a2_try) * sc.eta_bar);
            regime_ok = cond1 && cond2;
            if params.verbose && ~regime_ok
                fprintf('  Regime case1 fails: cond1=%d, cond2=%d\n', cond1, cond2);
                fprintf('    lam_bar=%.4f, (p/w_a2)*bar_eta=%.4f\n', ...
                        sc.lam_bar, (sc.p/w_a2_try)*sc.eta_bar);
            end

        case 'case2'
            cond1     = (hat_lam_fixed < sc.lam_bar);
            cond2     = (sc.lam_bar >= (sc.p / w_a2_try) * sc.eta_bar);
            regime_ok = cond1 && cond2;
            if params.verbose && ~regime_ok
                fprintf('  Regime case2 fails: cond1=%d, cond2=%d\n', cond1, cond2);
                fprintf('    lam_bar=%.4f, (p/w_a2)*bar_eta=%.4f\n', ...
                        sc.lam_bar, (sc.p/w_a2_try)*sc.eta_bar);
            end
    end

    if ~regime_ok
        if params.verbose
            fprintf('  Regime %s fails post-solution consistency check. Skipping.\n', regime);
        end
        continue
    end


    sol_found  = true;
    sol_regime = regime;
    v_sol      = v_sol_try;

    if params.verbose
        fprintf('  Regime %s ACCEPTED  (|R| = %.2e, w_a2 = %.6f)\n', ...
                regime, abs(R_check), w_a2_try);
    end
    break

end  % for regime



% no solution
if ~sol_found
    warning('solve_minwage_bgp: no regime converged for w_min = %.6f.', w_min);
    sol           = struct();
    sol.converged = false;
    sol.regime    = 'failed';
    sol.w_min     = w_min;
    return
end


w_a2 = exp(v_sol);

[Ia1, Ia2, Ix] = compute_integrals(w_min, w_a2, sol_regime, sc);
[Sigma, ~]     = compute_Sigma_Omega(Ia1, Ia2, Ix, w_min, w_a2, sc);

z_exp_y = sc.s_a2 * Sigma * (w_a2 ^ (1 + gam)) / (psi_A * Ia2);
z       = z_exp_y ^ (1 / sc.exp_y);

sol_tmp.z          = z;
sol_tmp.w_a1       = w_min;          % pinned at the floor
sol_tmp.w_a2       = w_a2;
sol_tmp.hat_eta    = w_a2 / w_min;
sol_tmp.hat_lambda = sc.p / w_min;   % the same threshold as above
sol_tmp.I_a1       = Ia1;
sol_tmp.I_a2       = Ia2;
sol_tmp.I_x        = Ix;
sol_tmp.Sigma      = Sigma;
sol_tmp.Omega      = sc.f * Sigma;
sol_tmp.converged  = true;
sol_tmp.regime_ok  = true;
sol_tmp.regime     = sol_regime;
sol = recover_post_solution(sol_tmp, sc, params);


L_a1 = sc.psi * sol.y * (w_min ^ (-1 - gam)) * Ia1 / Sigma;
u_a1 = sc.s_a1 - L_a1;

sol.L_a1  = L_a1;
sol.u_a1  = u_a1;
sol.w_min = w_min;


chk = sol.checks;

chk.u_a1_nonneg = (u_a1 >= 0);
if ~chk.u_a1_nonneg
    if params.verbose
        fprintf('FAILED: u_a1 = %.6f < 0. Min wage not binding. Returning failed.\n', u_a1);
    end
    sol.checks    = chk;
    sol.converged = false;
    sol.regime    = 'failed';
    return
end

chk.w_a2_above_floor = (w_a2 > w_min);
if ~chk.w_a2_above_floor
    if params.verbose
        fprintf('FAILED: w_a2 = %.6f <= w_min = %.6f. Joint-binding territory. Returning failed.\n', ...
                w_a2, w_min);
    end
    sol.checks    = chk;
    sol.converged = false;
    sol.regime    = 'failed';
    return
end

sol.checks    = chk;
sol.converged = true;

% rep
if params.verbose
    fprintf('\n--- Min-Wage BGP Solution ---\n');
    fprintf('  Regime       : %s (%s)\n', sol_regime, regime_labels.(sol_regime));
    fprintf('  w_min (=w_a1): %.6f\n', w_min);
    fprintf('  w_a2         : %.6f\n', w_a2);
    fprintf('  z            : %.6f\n', z);
    fprintf('  y            : %.6f\n', sol.y);
    fprintf('  k            : %.6f\n', sol.k);
    fprintf('  w_n          : %.6f\n', sol.w_n);
    fprintf('  X (robots)   : %.6f\n', sol.X);
    fprintf('  L_a1         : %.6f  (employed type-1)\n', L_a1);
    fprintf('  u_a1         : %.6f  (unemployed type-1)\n', u_a1);
    fprintf('  I_a1         : %.6f\n', Ia1);
    fprintf('  I_a2         : %.6f\n', Ia2);
    fprintf('  I_x          : %.6f\n', Ix);
    fprintf('\n--- Consistency Checks ---\n');
    fn = fieldnames(chk);
    for ii = 1:numel(fn)
        if chk.(fn{ii})
            status = 'PASS';
        else
            status = 'FAIL ***';
        end
        fprintf('  %-28s : %s\n', fn{ii}, status);
    end
    fprintf('\n');
end

end
