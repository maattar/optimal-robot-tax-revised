function sol = solve_minwage_bgp(params, w_min)
% SOLVE_MINWAGE_BGP  BGP with the wage floor binding for type-1 workers.
%
%   sol = solve_minwage_bgp(params, w_min)
%
%   Fixing w_a1 at the floor changes the problem in three ways.  First,
%   the automation margin becomes exogenous: hat_lambda = p/w_min is a
%   parameter comparison, settled before any solve, so only the choice
%   between the two automation cases is left to the numerics.  Second,
%   the system collapses to a single equation in w_a2, handled by fzero
%   rather than fsolve, with z recovered from (E3) afterwards.  Third,
%   type-1 market clearing is no longer imposed, so u_a1 = s_a1 - L_a1
%   comes out as an outcome rather than a condition.
%
%   This routine does not solve the competitive BGP to check the wage 
%   regime.
%
%   sol.regime is 'no_auto', 'case1', 'case2' or 'failed'.  It is 'failed'
%   whenever no regime converges, or the accepted solution implies
%   negative unemployment, or w_a2 falls to the floor.

sc    = compute_scalars(params);
psi_A = sc.psi * sc.A;
gam   = sc.gamma;

% --- Automation margin, settled exogenously ------------------------------
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

% Only the regimes compatible with that check are worth trying.
if auto_on
    regime_list = {'case1', 'case2'};
else
    regime_list = {'no_auto'};
end

regime_labels.no_auto = 'no-automation';
regime_labels.case1   = 'Case 1 (ray exits top)';
regime_labels.case2   = 'Case 2 (ray exits right)';

fzero_opts      = optimset('TolX', 1e-12, 'Display', 'off');
RESID_TOL       = 1e-8;
MAX_BRACK_ITERS = 40;

% w_a2 must exceed the floor, so the bracket scan starts just above it.
v_lo_base = log(w_min * 1.001);

% --- Try each admissible regime ------------------------------------------
sol_found  = false;
sol_regime = '';
v_sol      = NaN;

for ri = 1:numel(regime_list)
    regime = regime_list{ri};

    if params.verbose
        fprintf('\nTrying regime: %s ...\n', regime);
    end

    res_fun = @(v) residual_minwage(v, w_min, regime, sc);

    % Scan upward, doubling w_a2 each step, until two adjacent points give
    % opposite-sign residuals that are not both penalties.  Doubling covers
    % a wide range in few evaluations.
    bracket_ok = false;
    v_prev = v_lo_base;
    R_prev = res_fun(v_prev);

    for iter = 1:MAX_BRACK_ITERS
        v_next = v_prev + log(2);
        R_next = res_fun(v_next);

        if abs(R_prev) < 1e5 && abs(R_next) < 1e5 ...
                && sign(R_prev) ~= sign(R_next)
            v_lo       = v_prev;
            v_hi       = v_next;
            bracket_ok = true;
            break
        end

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

    % Check the residual directly as well: fzero can stop on a flat
    % stretch with an exit flag that looks fine.
    R_check = res_fun(v_sol_try);
    if abs(R_check) > RESID_TOL
        if params.verbose
            fprintf('  |R| = %.2e > %.2e at fzero solution. Skipping regime %s.\n', ...
                    abs(R_check), RESID_TOL, regime);
        end
        continue
    end

    % The solve assumed a regime; the solution has to satisfy it.
    w_a2_try = exp(v_sol_try);

    regime_ok = false;
    switch regime
        case 'no_auto'
            % The automation side already holds by the exogenous check, so
            % all that is left is that type-2 still occupies some tasks.
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

end

if ~sol_found
    warning('solve_minwage_bgp: no regime converged for w_min = %.6f.', w_min);
    sol           = struct();
    sol.converged = false;
    sol.regime    = 'failed';
    sol.w_min     = w_min;
    return
end

% --- Recover z from (E3) --------------------------------------------------
w_a2 = exp(v_sol);

[Ia1, Ia2, Ix] = compute_integrals(w_min, w_a2, sol_regime, sc);
[Sigma, ~]     = compute_Sigma_Omega(Ia1, Ia2, Ix, w_min, w_a2, sc);

z_exp_y = sc.s_a2 * Sigma * (w_a2 ^ (1 + gam)) / (psi_A * Ia2);
z       = z_exp_y ^ (1 / sc.exp_y);

% --- Remaining aggregates ------------------------------------------------
sol_tmp.z          = z;
sol_tmp.w_a1       = w_min;
sol_tmp.w_a2       = w_a2;
sol_tmp.hat_eta    = w_a2 / w_min;
sol_tmp.hat_lambda = sc.p / w_min;
sol_tmp.I_a1       = Ia1;
sol_tmp.I_a2       = Ia2;
sol_tmp.I_x        = Ix;
sol_tmp.Sigma      = Sigma;
sol_tmp.Omega      = sc.f * Sigma;
sol_tmp.converged  = true;
sol_tmp.regime_ok  = true;
sol_tmp.regime     = sol_regime;
sol = recover_post_solution(sol_tmp, sc, params);

% --- Type-1 employment and unemployment ----------------------------------
% Labour demand at the floor, against a supply of s_a1 that is no longer
% required to clear.  In the competitive BGP the two coincide.
L_a1 = sc.psi * sol.y * (w_min ^ (-1 - gam)) * Ia1 / Sigma;
u_a1 = sc.s_a1 - L_a1;

sol.L_a1  = L_a1;
sol.u_a1  = u_a1;
sol.w_min = w_min;

% --- Checks specific to the floor ----------------------------------------
chk = sol.checks;

% Negative unemployment means the floor never bound in the first place and
% the competitive allocation still clears; return failure so the caller
% penalises the draw rather than accepting a degenerate solution.
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

% w_a2 at or below the floor is the joint-binding case, outside this
% variant.  The bracket starts above w_min so it should not arise, but the
% check is kept explicit so no caller can receive a converged solution
% from an inadmissible regime.
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
