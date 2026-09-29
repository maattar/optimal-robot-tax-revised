function R = residual_minwage(log_wa2, w_min, regime, sc)

%   R = residual_minwage(log_wa2, w_min, regime, sc)

w_a2 = exp(log_wa2);

[Ia1, Ia2, Ix] = compute_integrals(w_min, w_a2, regime, sc);

% regime check
if Ia2 <= 0 || Ix < 0
    R = 1e6;
    return
end


[Sigma, Omega] = compute_Sigma_Omega(Ia1, Ia2, Ix, w_min, w_a2, sc);

psi_A     = sc.psi * sc.A;    
s_a2      = sc.s_a2;          % the supply of type-2 automatable labour
gam       = sc.gamma;         % theta/(1-theta), between -1 and 0
exp_psi   = sc.exp_y;         % psi/(1-alpha)
exp_z     = sc.exp_z;         % (1-alpha-psi)/(1-alpha)
exp_Omega = sc.exp_Omega;     % (1-theta)/theta, negative

kappa = exp_z / exp_psi;      % (1-alpha-psi)/psi  > 0


inner_LHS = s_a2 * Sigma * (w_a2 ^ (1 + gam)) / (psi_A * Ia2);
LHS       = inner_LHS ^ kappa;


RHS = psi_A * (Omega ^ exp_Omega);


R = LHS / RHS - 1;

end
