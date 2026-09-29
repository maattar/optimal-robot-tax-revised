function R = residual_minwage(log_wa2, w_min, regime, sc)
% RESIDUAL_MINWAGE  Residual of the collapsed min-wage BGP equation.
%
%   R = residual_minwage(log_wa2, w_min, regime, sc)
%
%   With w_a1 pinned at the floor, (E1) and (E3) can be reduced to one
%   equation in w_a2.  Solving (E3) for z^exp_y and raising it to
%   kappa = exp_z/exp_y = (1-alpha-psi)/psi gives z^exp_z, which
%   substituted into (E1) leaves
%
%     [s_a2*Sigma*w_a2^(1+gamma) / (psi_A*Ia2)]^kappa = psi_A*(f*Sigma)^exp_Omega.
%
%   The residual is the ratio of the two sides minus one.  w_a2 is passed
%   in logs to keep the search on the positive half-line, and a regime
%   mismatch returns a large penalty so that fzero stays clear of it.

w_a2 = exp(log_wa2);

[Ia1, Ia2, Ix] = compute_integrals(w_min, w_a2, regime, sc);

% Ia2 <= 0 means type-2 has been priced out entirely (hat_eta has reached
% eta_bar); Ix < 0 means the robot region has collapsed.  Either way the
% guess has left the assumed regime.
if Ia2 <= 0 || Ix < 0
    R = 1e6;
    return
end

[Sigma, Omega] = compute_Sigma_Omega(Ia1, Ia2, Ix, w_min, w_a2, sc);

psi_A     = sc.psi * sc.A;
s_a2      = sc.s_a2;
gam       = sc.gamma;
exp_psi   = sc.exp_y;
exp_z     = sc.exp_z;
exp_Omega = sc.exp_Omega;

kappa = exp_z / exp_psi;

inner_LHS = s_a2 * Sigma * (w_a2 ^ (1 + gam)) / (psi_A * Ia2);
LHS       = inner_LHS ^ kappa;

RHS = psi_A * (Omega ^ exp_Omega);

R = LHS / RHS - 1;

end
