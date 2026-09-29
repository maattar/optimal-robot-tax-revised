function R = residuals_3D(v, sc, regime)
% RESIDUALS_3D  Residuals of the three-equation competitive BGP system.
%
%   R = residuals_3D(v, sc, regime)
%
%   The unknowns are carried in logs, v = [log z; log w_a1; log w_a2], so
%   that fsolve cannot step into negative prices and no explicit bounds
%   are needed.  regime is 'case1' or 'case2'.
%
%   The equations are (E1) and the two labour-market clearing conditions
%   for the automatable types; each residual is written as a ratio minus
%   one, which makes all three dimensionless and of comparable size near
%   the solution.  A regime mismatch shows up as a negative integral and
%   returns a large penalty.

    z    = exp(v(1));
    w_a1 = exp(v(2));
    w_a2 = exp(v(3));

    [Ia1, Ia2, Ix] = compute_integrals(w_a1, w_a2, regime, sc);

    PENALTY = 1e6;
    if Ia2 < 0 || Ix < 0
        R = PENALTY * ones(3, 1);
        return
    end

    [Sigma, Omega] = compute_Sigma_Omega(Ia1, Ia2, Ix, w_a1, w_a2, sc);

    if Sigma <= 0 || Omega <= 0
        R = PENALTY * ones(3, 1);
        return
    end

    psi_A  = sc.psi * sc.A;
    exp_y  = sc.exp_y;
    exp_z  = sc.exp_z;
    exp_Om = sc.exp_Omega;
    gamma  = sc.gamma;
    s_a1   = sc.s_a1;
    s_a2   = sc.s_a2;

    % (E1) conditional demand for the intermediate.
    lhs_E1 = z^exp_z;
    rhs_E1 = psi_A * Omega^exp_Om;
    R1     = lhs_E1 / rhs_E1 - 1;

    % (E2) type-1 labour market.
    lhs_E2 = psi_A * z^exp_y * w_a1^(-1 - gamma) * Ia1 / Sigma;
    R2     = lhs_E2 / s_a1 - 1;

    % (E3) type-2 labour market.
    lhs_E3 = psi_A * z^exp_y * w_a2^(-1 - gamma) * Ia2 / Sigma;
    R3     = lhs_E3 / s_a2 - 1;

    R = [R1; R2; R3];

end
