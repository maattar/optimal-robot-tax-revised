function [Ia1, Ia2, Ix] = compute_integrals(w_a1, w_a2, regime, sc)
% COMPUTE_INTEGRALS  Task-space integrals for the three factor types.
%
%   [Ia1, Ia2, Ix] = compute_integrals(w_a1, w_a2, regime, sc)
%
%   Integrates the density-weighted task shares over the rectangle
%   (1, eta_bar] x (0, lam_bar].  Which regions belong to which factor is
%   fixed by two thresholds,
%       lam_hat = p / w_a1      robots undercut type-1 above this lambda
%       eta_hat = w_a2 / w_a1   type-2 undercuts type-1 above this eta
%   and by the dividing ray lam = (p/w_a2)*eta between type-2 and robots.
%   regime is 'no_auto', 'case1' or 'case2'.
%
%   The integrals are non-negative when the assumed regime is the right
%   one.  A negative Ia2 or Ix during a solve means the iterate has
%   crossed a regime boundary, and the callers use that as a signal.

p       = sc.p;
gam     = sc.gamma;   %#ok<NASGU>
g1      = sc.g1;      % gamma + 1
g2      = sc.g2;      % gamma + 2
eta_bar = sc.eta_bar;
lam_bar = sc.lam_bar;

lam_hat = p   / w_a1;
eta_hat = w_a2 / w_a1;

switch regime

    case 'no_auto'
    % Robots priced out (lam_hat >= lam_bar).  The whole rectangle is
    % split vertically at eta_hat, type-1 to the left, type-2 to the right.

        Ia1 = lam_bar * (eta_hat - 1);

        Ia2 = lam_bar * (eta_bar^g1 - eta_hat^g1) / g1;

        Ix  = 0;

    case 'case1'
    % The ray leaves the rectangle through the top edge, at eta_bar'
    % below eta_bar.  Type-2 therefore splits in two: under the ray up to
    % eta_prime, then a full-height strip out to eta_bar.  Type-1 keeps
    % the rectangle [1, eta_hat] x [0, lam_hat].

        eta_prime = (w_a2 / p) * lam_bar;

        Ia1 = lam_hat * (eta_hat - 1);

        Ia2 = (p / w_a2) * (eta_prime^g2 - eta_hat^g2) / g2 ...
            + lam_bar    * (eta_bar^g1    - eta_prime^g1) / g1;

        Ix  = (w_a2 / p) * (lam_bar^g2 - lam_hat^g2) / g2 ...
            -              (lam_bar^g1 - lam_hat^g1) / g1;

    case 'case2'
    % The ray leaves through the right edge instead, at lam_prime below
    % lam_bar.  Now it is the robot region that splits: bounded by the ray
    % up to lam_prime, then the full width of the rectangle above it.
    % Type-2 needs no split.

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
