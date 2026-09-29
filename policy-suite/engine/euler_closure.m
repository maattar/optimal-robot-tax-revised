function r = euler_closure(Phi, tau_K)
% EULER_CLOSURE  The interest rate the household Euler equation implies at
% a given capital tax rate.
%
%   r = euler_closure(Phi, tau_K)
%
% At tau_K = 1 no BGP exists and r runs off to infinity below it.

TAU_K_SING_EPS = 1e-8;

if ~isscalar(Phi) || ~isfinite(Phi) || ~isreal(Phi) || Phi <= 0
    error(['euler_closure: Phi must be a finite positive real scalar ' ...
           '(got a value that is not). Phi = rho + g_N + sigma*g_A is ' ...
           'the household''s after-tax required return and is positive ' ...
           'in any sensible calibration; a non-positive Phi means the ' ...
           'loaded calibration .mat is not the one this suite expects.']);
end

r = Phi ./ (1 - tau_K);

r(tau_K >= 1 - TAU_K_SING_EPS) = Inf;

end
