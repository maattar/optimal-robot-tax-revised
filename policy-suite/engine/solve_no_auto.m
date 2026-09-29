function sol = solve_no_auto(sc, params)
%SOLVE_NO_AUTO  blanced growth where nothing is automated.
%
%   sol = solve_no_auto(sc, params)
%
% one equation in hat_eta:
%
%   s_a1/s_a2 = hat_eta^(1+gamma) * (gamma+1)*(hat_eta-1)
%               / (eta_bar^(gamma+1) - hat_eta^(gamma+1))
%

    gamma   = sc.gamma;          % theta/(1-theta) in (-1,0)
    g1      = sc.g1;             
    f       = sc.f;
    psi_A   = sc.psi * sc.A;    % psi*A
    exp_y   = sc.exp_y;          % psi/(1-alpha)
    exp_z   = sc.exp_z;          % (1-alpha-psi)/(1-alpha)
    exp_Om  = sc.exp_Omega;      % (1-theta)/theta [< 0]c
    p       = sc.p;

    eta_bar = sc.eta_bar;
    lam_bar = sc.lam_bar;        
    s_a1    = sc.s_a1;
    s_a2    = sc.s_a2;

    

    ratio = s_a1 / s_a2;

    % The residual
    resid_NA = @(eta) ...
        eta^(1+gamma) * g1 * (eta - 1) / (eta_bar^g1 - eta^g1) - ratio;

    
    eta_lo = 1 + 1e-10;
    eta_hi = eta_bar * (1 - 1e-10);

    val_lo = resid_NA(eta_lo);
    val_hi = resid_NA(eta_hi);

    if val_lo * val_hi >= 0
        warning('solve_no_auto: fzero bracket has no sign change. Check parameters (eta_bar=%.4g, s_a1/s_a2=%.4g).', eta_bar, ratio);
        sol = make_empty_sol();
        return
    end

    opts    = optimset('TolX', 1e-12, 'Display', 'off');
    hat_eta = fzero(resid_NA, [eta_lo, eta_hi], opts);

    % the integrals

    Ia1 = lam_bar * (hat_eta - 1);
    Ia2 = lam_bar * (eta_bar^g1 - hat_eta^g1) / g1;
    Ix  = 0;

    
    B     = Ia1 + hat_eta^(-gamma) * Ia2;

    exp_r = exp_z / exp_y;      % (1-alpha-psi)/psi, positive
    exp_w = exp_r + 1;          % (1-alpha)/psi, positive

    K    = psi_A ...
           * (psi_A * Ia1 / s_a1)^exp_r ...
           * f^exp_Om ...
           * B^(exp_Om - exp_r);

    w_a1 = K^(1/exp_w);        % = K^{psi/(1-alpha)} = K^{exp_y}

    %
    % =========================================================================

    w_a2 = hat_eta * w_a1;

   
    z = (s_a1 * B * w_a1 / (psi_A * Ia1))^(1/exp_y);

  

  [Sigma, Omega] = compute_Sigma_Omega(Ia1, Ia2, Ix, w_a1, w_a2, sc);

    hat_lambda = p / w_a1;

    cond_no_auto = (hat_lambda >= lam_bar);
    cond_eta_hi  = (hat_eta < eta_bar);
    cond_eta_lo  = (hat_eta > 1);

    regime_ok = cond_no_auto && cond_eta_hi && cond_eta_lo;

    if ~cond_no_auto
        warning('solve_no_auto: regime check FAILED (hat_lambda=%.4g < lam_bar=%.4g). Solution likely belongs to an automation regime.', hat_lambda, lam_bar);
    end
    if ~cond_eta_hi
        warning('solve_no_auto: hat_eta=%.4g >= eta_bar=%.4g. Type-2 workers priced out of all tasks.', hat_eta, eta_bar);
    end

    %output
    % =========================================================================

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
    sol.regime     = 'no_auto';

end

% =========================================================================
function sol = make_empty_sol()
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
    sol.regime     = 'no_auto';
end
