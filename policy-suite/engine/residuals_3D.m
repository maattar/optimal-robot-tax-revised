function R = residuals_3D(v, sc, regime)
%RESIDUALS_3D  The residuals of the three-equation competitive system
%
%   R = residuals_3D(v, sc, regime)
%
% The unknowns are v = [log z; log w_a1; log w_a2]. 

    z    = exp(v(1));
    w_a1 = exp(v(2));
    w_a2 = exp(v(3));

    % The task-space integrals
    [Ia1, Ia2, Ix] = compute_integrals(w_a1, w_a2, regime, sc);

    % regime check
    PENALTY = 1e6;
    if Ia2 < 0 || Ix < 0
        R = PENALTY * ones(3, 1);
        return
    end

    % The cost aggregate
    [Sigma, Omega] = compute_Sigma_Omega(Ia1, Ia2, Ix, w_a1, w_a2, sc);

    if Sigma <= 0 || Omega <= 0
        R = PENALTY * ones(3, 1);
        return
    end

    psi_A  = sc.psi * sc.A;
    exp_y  = sc.exp_y;         % psi/(1-alpha)
    exp_z  = sc.exp_z;         % (1-alpha-psi)/(1-alpha)
    exp_Om = sc.exp_Omega;     % (1-theta)/theta, negative
    gamma  = sc.gamma;         % theta/(1-theta), between -1 and 0
    s_a1   = sc.s_a1;
    s_a2   = sc.s_a2;

    % Conditional factor demand
    lhs_E1 = z^exp_z;
    rhs_E1 = psi_A * Omega^exp_Om;
    R1     = lhs_E1 / rhs_E1 - 1;

    % Clearing in the type-1 labour market
    lhs_E2 = psi_A * z^exp_y * w_a1^(-1 - gamma) * Ia1 / Sigma;
    R2     = lhs_E2 / s_a1 - 1;

    % Clearing in the type-2 labour market
    lhs_E3 = psi_A * z^exp_y * w_a2^(-1 - gamma) * Ia2 / Sigma;
    R3     = lhs_E3 / s_a2 - 1;

    R = [R1; R2; R3];

end
