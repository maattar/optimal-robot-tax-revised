function [tauR_c, branch, w_a1_eff] = corner_of_pol(r, delta, lam_bar, w_min, C, alpha)
% CORNER_OF_POL  The automation-blocking corner, with the wage floor taken
% as an argument rather than read from the calibration. One formula serves
% every policy.
%
%   [tauR_c, branch, w_a1_eff] = corner_of_pol(r, delta, lam_bar, w_min, C, alpha)
%
%       tau_R^corner = lam_bar * w_a1_eff / (r+delta) - 1
%       w_a1_eff     = max( w_min , w_a1^{c,NA} )
%       w_a1^{c,NA}  = C * (r+delta)^(-alpha/(1-alpha))
%
% The floor is an argument because reading it from a cp struct is right
% only while it stays pinned for a whole search, as in pol4. pol6 frees
% it, so the corner has to be evaluated at each iterate's own floor, and
% passing it in removes the temptation: there is no cp.w_min to reach
% for.



if ~isfinite(r)
    tauR_c = NaN; branch = 'singular'; w_a1_eff = NaN;
    return
end

rd       = r + delta;
w_a1_cNA = C * rd^(-alpha / (1 - alpha));

if w_a1_cNA >= w_min
    w_a1_eff = w_a1_cNA;
    branch   = 'slack';
else
    w_a1_eff = w_min;
    branch   = 'binding';
end

tauR_c = lam_bar * w_a1_eff / rd - 1;

end
