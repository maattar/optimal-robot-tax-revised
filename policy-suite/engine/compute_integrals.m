function [Ia1, Ia2, Ix] = compute_integrals(w_a1, w_a2, regime, sc)
% COMPUTE_INTEGRALS  The three task-space integrals the BGP
% solver needs, each the area-weighted share of tasks going to one factor.
%
%   [Ia1, Ia2, Ix] = compute_integrals(w_a1, w_a2, regime, sc)

p       = sc.p;
gam     = sc.gamma;   %#ok<NASGU>  not used directly; g1 and g2 are
g1      = sc.g1;      % gamma + 1
g2      = sc.g2;      % gamma + 2
eta_bar = sc.eta_bar;
lam_bar = sc.lam_bar;

% The two thresholds that allocate tasks
lam_hat = p   / w_a1;   % robots take a task only above this in lambda
eta_hat = w_a2 / w_a1;  % type-2 labour takes one only above this in eta

% Take the regime as given and follow it
switch regime

    case 'no_auto'
    
        Ia1 = lam_bar * (eta_hat - 1);

        Ia2 = lam_bar * (eta_bar^g1 - eta_hat^g1) / g1;

        Ix  = 0;

    case 'case1'

        eta_prime = (w_a2 / p) * lam_bar;

        Ia1 = lam_hat * (eta_hat - 1);

        Ia2 = (p / w_a2) * (eta_prime^g2 - eta_hat^g2) / g2 ...
            + lam_bar    * (eta_bar^g1    - eta_prime^g1) / g1;

        Ix  = (w_a2 / p) * (lam_bar^g2 - lam_hat^g2) / g2 ...
            -              (lam_bar^g1 - lam_hat^g1) / g1;

    case 'case2'

        lam_prime = (p / w_a2) * eta_bar;

        Ia1 = lam_hat * (eta_hat - 1);

        Ia2 = (p / w_a2) * (eta_bar^g2 - eta_hat^g2) / g2;

        Ix  = (w_a2 / p) * (lam_prime^g2 - lam_hat^g2) / g2 ...
            -              (lam_prime^g1 - lam_hat^g1) / g1 ...
            + (eta_bar - 1) / g1 * (lam_bar^g1 - lam_prime^g1);

    otherwise
        error('compute_integrals: unknown regime "%s". Valid: no_auto | case1 | case2.', ...
              regime);

end

end
