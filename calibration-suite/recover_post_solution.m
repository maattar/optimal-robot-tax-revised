function out = recover_post_solution(sol, sc, params)
% RECOVER_POST_SOLUTION  Back out y, k, w_n and X, and check the solution.
%
%   out = recover_post_solution(sol, sc, params)
%
%   Once z and the automatable wages are known, the remaining aggregates
%   follow directly: y from the reduced form, k from the capital first-
%   order condition, w_n from the residual claim on output, and robot
%   demand X from its share of the intermediate bill.  X is zero by
%   construction under no automation.
%
%   The returned struct is sol with y, k, w_n, X and a checks field
%   appended.  Failed checks raise a warning rather than an error.

    z      = sol.z;
    w_a1   = sol.w_a1;
    w_a2   = sol.w_a2;
    Sigma  = sol.Sigma;
    I_x    = sol.I_x;
    regime = sol.regime;

    A       = sc.A;
    exp_y   = sc.exp_y;
    gamma   = sc.gamma;
    p       = sc.p;
    alpha   = sc.alpha;
    psi     = sc.psi;
    r       = sc.r;
    delta   = sc.delta;
    s_n     = sc.s_n;
    eta_bar = sc.eta_bar;
    lam_bar = sc.lam_bar;

    verbose = isfield(params, 'verbose') && params.verbose;

    % --- Recovery ---------------------------------------------------------

    y   = A * z^exp_y;
    k   = (alpha / (r + delta)) * y;
    w_n = (1 - alpha - psi) * y / s_n;

    if strcmp(regime, 'no_auto')
        X = 0;
    else
        X = (psi * y / p) * (p^(-gamma) * I_x / Sigma);
    end

    % --- Checks -----------------------------------------------------------

    checks = struct();

    checks.eta_gt_1   = sol.hat_eta > 1;
    checks.eta_lt_bar = sol.hat_eta < eta_bar;

    if strcmp(regime, 'no_auto')
        checks.lambda_ok = true;
    else
        checks.lambda_ok = sol.hat_lambda < lam_bar;
    end

    checks.wages_pos = (w_a1 > 0) && (w_a2 > 0) && (w_n > 0);
    checks.X_nonneg  = (X >= 0);

    checks.all_ok = checks.eta_gt_1 && checks.eta_lt_bar && ...
                    checks.lambda_ok && checks.wages_pos && checks.X_nonneg;

    if verbose
        fprintf('\n  [POST-SOL] Consistency checks:\n');
        fprintf('    eta_gt_1   (hat_eta > 1):        %s  (hat_eta=%.4g)\n', ...
                yn(checks.eta_gt_1),   sol.hat_eta);
        fprintf('    eta_lt_bar (hat_eta < eta_bar):  %s  (eta_bar=%.4g)\n', ...
                yn(checks.eta_lt_bar), eta_bar);
        if ~strcmp(regime, 'no_auto')
            fprintf('    lambda_ok  (hat_lam < lam_bar): %s  (hat_lambda=%.4g, lam_bar=%.4g)\n', ...
                    yn(checks.lambda_ok), sol.hat_lambda, lam_bar);
        end
        fprintf('    wages_pos:                        %s  (w_a1=%.4g  w_a2=%.4g  w_n=%.4g)\n', ...
                yn(checks.wages_pos), w_a1, w_a2, w_n);
        fprintf('    X_nonneg:                         %s  (X=%.4g)\n', ...
                yn(checks.X_nonneg), X);
        fprintf('    --- all_ok: %s ---\n', yn(checks.all_ok));

        wages  = [w_a1, w_a2, w_n];
        labels = {'w_a1', 'w_a2', 'w_n'};
        [~, ord] = sort(wages);
        fprintf('    Wage ordering: ');
        for i = 1:3
            fprintf('%s=%.4g', labels{ord(i)}, wages(ord(i)));
            if i < 3, fprintf(' < '); end
        end
        fprintf('\n');
    end

    if ~checks.all_ok
        warning('recover_post_solution:checkFailed', ...
            'One or more post-solution consistency checks failed. Inspect out.checks.');
    end

    out        = sol;
    out.y      = y;
    out.k      = k;
    out.w_n    = w_n;
    out.X      = X;
    out.checks = checks;

end

function s = yn(flag)
    if flag, s = 'OK  '; else, s = 'FAIL'; end
end
