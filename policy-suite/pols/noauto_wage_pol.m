function [wa1, C] = noauto_wage_pol(r, cp)
% NOAUTO_WAGE_POL  The competitive type-1 wage when nothing is automated,
% at a given interest rate, with the scale constant C in
%
%     w_a1^{c,NA} = C * (r+delta)^(-alpha/(1-alpha))
%
%   [wa1, C] = noauto_wage_pol(r, cp)
%
% C does not vary with r. The rate is taken directly rather than derived
% from tau_K, so one file serves pol6, where r is fixed, and pol4, where
% the caller has resolved r through euler_closure first.
%

%
%   cp   alpha, psi, theta, lam_bar, eta_bar, delta, s_a1, s_a2.
p = cp;
p.tau_R   = 0;
p.w_min   = 0;
p.r       = r;
p.s_n     = 1 - p.s_a1 - p.s_a2;
p.verbose = false;

sc = compute_scalars(p);
ws = warning('off', 'all');    
sol = solve_no_auto(sc, p);
warning(ws);

wa1 = sol.w_a1;
C   = wa1 * (r + cp.delta)^(cp.alpha / (1 - cp.alpha));

end
