function [Sigma, Omega] = compute_Sigma_Omega(I_a1, I_a2, I_x, w_a1, w_a2, sc)
% COMPUTE_SIGMA_OMEGA  Cost aggregate of the task technology.
%
%   Sigma weights the three task integrals by their factor prices and is
%   the only place where factor prices enter the reduced form: it sets the
%   scale of z through (E1), and each labour-market clearing condition
%   carries the share term w_ai^(-1-gamma)*I_ai/Sigma.  Omega = f*Sigma is
%   the form in which it appears in (E1).
%
%   Since gamma is negative, every w^(-gamma) is positive and Sigma is
%   positive whenever one integral is.  A non-positive Sigma means either
%   a wrong regime or an infeasible wage guess.

    gamma = sc.gamma;
    p     = sc.p;
    f     = sc.f;

    % The robot term drops out on its own under no automation, where I_x = 0.
    Sigma = w_a1^(-gamma) * I_a1 ...
          + w_a2^(-gamma) * I_a2 ...
          + p^(-gamma)    * I_x;

    Omega = f * Sigma;

end
