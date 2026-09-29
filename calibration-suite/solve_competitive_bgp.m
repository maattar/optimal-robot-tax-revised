function sol = solve_competitive_bgp(sc, params)
% SOLVE_COMPETITIVE_BGP  Competitive BGP, dispatching over regimes.
%
%   sol = solve_competitive_bgp(sc, params)
%
%   The regime is not known in advance, so the three candidates are tried
%   in order — no automation, then the two automation cases — and the
%   first one that both converges and reproduces its own defining
%   inequalities is returned.  Convergence means a positive fsolve flag
%   and a residual below RESID_TOL; the flag alone is not enough.
%
%   sol.regime is 'no_auto', 'case1', 'case2' or 'failed'.
%   Set params.verbose to trace the dispatch.

    RESID_TOL = 1e-8;
    verbose   = isfield(params, 'verbose') && params.verbose;

    % --- No automation ---------------------------------------------------

    sol = solve_no_auto(sc, params);

    if sol.converged && sol.regime_ok
        if verbose
            fprintf('  [BGP] Regime confirmed: no-automation.\n');
            fprintf('        hat_lambda=%.4g >= lam_bar=%.4g\n', ...
                    sol.hat_lambda, sc.lam_bar);
        end
        return
    end

    if verbose && sol.converged
        fprintf('  [BGP] No-auto solve converged but regime check failed\n');
        fprintf('        (hat_lambda=%.4g < lam_bar=%.4g). Trying automation.\n', ...
                sol.hat_lambda, sc.lam_bar);
    end

    % The no-automation wages are the natural warm start for the automation
    % cases; they sit just the other side of the boundary.  Failing that,
    % put w_a1 slightly under p/lam_bar so automation is at least active.
    if sol.converged
        v0 = [log(sol.z); log(sol.w_a1); log(sol.w_a2)];
    else
        wa1_fb = sc.p / sc.lam_bar * 0.8;
        wa2_fb = wa1_fb * max(1.5, 0.5 * (1 + sc.eta_bar));
        v0     = [log(1.0); log(wa1_fb); log(wa2_fb)];
        if verbose
            fprintf('  [BGP] No-auto solve failed; using fallback initial guess.\n');
        end
    end

    opts = optimoptions('fsolve', ...
        'TolFun',                1e-10, ...
        'TolX',                  1e-10, ...
        'MaxIterations',         1000,  ...
        'MaxFunctionEvaluations',5000,  ...
        'Display',               'off');

    % --- Case 1: dividing ray exits through the top ----------------------

    [v1, R1_vec, flag1] = fsolve(@(v) residuals_3D(v, sc, 'case1'), v0, opts);

    converged_c1 = (flag1 > 0) && (norm(R1_vec, Inf) < RESID_TOL);

    if converged_c1
        sol_c1 = unpack_auto_sol(v1, sc, 'case1');
        if verbose
            print_sol_summary(sol_c1, 'case1', R1_vec, sc);
        end
        if sol_c1.regime_ok
            sol = sol_c1;
            if verbose, fprintf('  [BGP] Regime confirmed: automation Case 1.\n'); end
            return
        else
            if verbose, fprintf('  [BGP] Case 1 converged but regime check failed. Trying Case 2.\n'); end
        end
        v0_c2 = v1;
    else
        if verbose
            fprintf('  [BGP] Case 1 fsolve did not converge (flag=%d). Trying Case 2.\n', flag1);
        end
        v0_c2 = v0;
    end

    % --- Case 2: dividing ray exits through the right --------------------

    [v2, R2_vec, flag2] = fsolve(@(v) residuals_3D(v, sc, 'case2'), v0_c2, opts);

    converged_c2 = (flag2 > 0) && (norm(R2_vec, Inf) < RESID_TOL);

    if converged_c2
        sol_c2 = unpack_auto_sol(v2, sc, 'case2');
        if verbose
            print_sol_summary(sol_c2, 'case2', R2_vec, sc);
        end
        if sol_c2.regime_ok
            sol = sol_c2;
            if verbose, fprintf('  [BGP] Regime confirmed: automation Case 2.\n'); end
            return
        else
            if verbose, fprintf('  [BGP] Case 2 converged but regime check failed.\n'); end
        end
    else
        if verbose
            fprintf('  [BGP] Case 2 fsolve did not converge (flag=%d).\n', flag2);
        end
    end

    warning('solve_competitive_bgp:noSolution', ...
        ['No self-consistent regime found.\n' ...
         'Check parameters, or try adjusting initial guesses.']);
    sol = make_failure_sol();

end


function sol = unpack_auto_sol(v, sc, regime)
% Rebuild the solution from the log-unknowns and test the inequalities
% that define the regime that was assumed during the solve.

    z    = exp(v(1));
    w_a1 = exp(v(2));
    w_a2 = exp(v(3));

    [Ia1, Ia2, Ix] = compute_integrals(w_a1, w_a2, regime, sc);
    [Sigma, Omega]  = compute_Sigma_Omega(Ia1, Ia2, Ix, w_a1, w_a2, sc);

    hat_eta    = w_a2 / w_a1;
    hat_lambda = sc.p / w_a1;
    lam_bar    = sc.lam_bar;
    eta_bar    = sc.eta_bar;
    p          = sc.p;

    % Automation must be active in both cases; the cases differ only in
    % which edge the dividing ray leaves through.
    auto_ok = hat_lambda < lam_bar;

    switch regime
        case 'case1'
            geom_ok = lam_bar < (p / w_a2) * eta_bar;
        case 'case2'
            geom_ok = lam_bar >= (p / w_a2) * eta_bar;
        otherwise
            geom_ok = false;
    end

    eta_ok    = hat_eta < eta_bar;
    regime_ok = auto_ok && geom_ok && eta_ok;

    sol.z          = z;
    sol.w_a1       = w_a1;
    sol.w_a2       = w_a2;
    sol.hat_eta    = hat_eta;
    sol.hat_lambda = hat_lambda;
    sol.I_a1       = Ia1;
    sol.I_a2       = Ia2;
    sol.I_x        = Ix;
    sol.Sigma      = Sigma;
    sol.Omega      = Omega;
    sol.converged  = true;
    sol.regime_ok  = regime_ok;
    sol.regime     = regime;
end


function sol = make_failure_sol()
    sol.z          = NaN;
    sol.w_a1       = NaN;
    sol.w_a2       = NaN;
    sol.hat_eta    = NaN;
    sol.hat_lambda = NaN;
    sol.I_a1       = NaN;
    sol.I_a2       = NaN;
    sol.I_x        = NaN;
    sol.Sigma      = NaN;
    sol.Omega      = NaN;
    sol.converged  = false;
    sol.regime_ok  = false;
    sol.regime     = 'failed';
end


function print_sol_summary(sol, regime_label, R_vec, sc)
    fprintf('  [BGP] %s solve: ||R||_inf = %.2e\n', regime_label, norm(R_vec, Inf));
    fprintf('        z=%.4g  w_a1=%.4g  w_a2=%.4g\n', sol.z, sol.w_a1, sol.w_a2);
    fprintf('        hat_lambda=%.4g (lam_bar=%.4g)  hat_eta=%.4g (eta_bar=%.4g)\n', ...
            sol.hat_lambda, sc.lam_bar, sol.hat_eta, sc.eta_bar);
    geom_val = (sc.p / sol.w_a2) * sc.eta_bar;
    fprintf('        (p/w_a2)*eta_bar=%.4g  lam_bar=%.4g\n', geom_val, sc.lam_bar);
end
