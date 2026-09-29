function [Sigma, Omega] = compute_Sigma_Omega(I_a1, I_a2, I_x, w_a1, w_a2, sc)
%COMPUTE_SIGMA_OMEGA  The task-weighted cost aggregate Sigma, and Omega.
%
%   [Sigma, Omega] = compute_Sigma_Omega(I_a1, I_a2, I_x, w_a1, w_a2, sc)
%
%   Sigma = w_a1^-gamma I_a1 + w_a2^-gamma I_a2 + p^-gamma I_x,  Omega = f*Sigma


    gamma = sc.gamma;   % between -1 and 0
    p     = sc.p;
    f     = sc.f;

    Sigma = w_a1^(-gamma) * I_a1 ...
          + w_a2^(-gamma) * I_a2 ...
          + p^(-gamma)    * I_x;

    Omega = f * Sigma;

end
