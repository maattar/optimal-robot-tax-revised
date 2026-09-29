function [Phi, rel_gap] = phi_from_calibration(fp, verbose)
% PHI_FROM_CALIBRATION  The Euler constant Phi out of a loaded calibration,
% re-checking consistency on every read.
%
%   [Phi, rel_gap] = phi_from_calibration(fp)
%   [Phi, rel_gap] = phi_from_calibration(fp, verbose)

if nargin < 2 || isempty(verbose)
    verbose = true;
end

need = {'r', 'tau_K', 'rho', 'g_N', 'g_A', 'sigma'};
missing = need(~cellfun(@(fn) isfield(fp, fn), need));
if ~isempty(missing)
    error(['phi_from_calibration: fp is missing required field(s): %s. Phi needs all ' ...
           'of r, tau_K, rho, g_N, g_A, sigma.'], strjoin(missing, ', '));
end

if fp.tau_K >= 1
    error(['phi_from_calibration: calibrated tau_K = %.6f >= 1, so r = Phi/(1-tau_K) ' ...
           'is undefined at the baseline itself.'], fp.tau_K);
end
if fp.g_A == 0
    error(['phi_from_calibration: g_A = 0, so sigma is not identified by the Euler ' ...
           'equation and the consistency check is undefined.']);
end

Phi_implied = fp.r * (1 - fp.tau_K);                 % the form to use
Phi_struct  = fp.rho + fp.g_N + fp.sigma * fp.g_A;   % for the check only

if Phi_implied <= 0
    error(['phi_from_calibration: Phi = fp.r*(1-fp.tau_K) = %.6g is non-positive. The ' ...
           'household''s after-tax required return must be positive.'], Phi_implied);
end

rel_gap = abs(Phi_struct - Phi_implied) / abs(Phi_implied);

TOL_EXACT = 1e-12;
TOL_ROUND = 1e-5;

if rel_gap < TOL_EXACT
    if verbose
        fprintf(['  Phi check: [PASS - EXACT]  rel gap %.2e.  ' ...
                 'Phi = %.11f\n'], rel_gap, Phi_implied);
    end
elseif rel_gap < TOL_ROUND
    warning(['phi_from_calibration: [PASS - WITH CAVEAT] Phi formulas differ by rel %.2e. ' ...
             'This is the signature of a ROUNDED stored sigma, not a genuine ' ...
             'inconsistency. Proceeding with Phi = fp.r*(1-fp.tau_K) = %.11f, ' ...
             'which is the correct choice precisely in this case: it ' ...
             'returns fp.r at tau_K_base and so preserves nesting.'], ...
            rel_gap, Phi_implied);
else
    error(['phi_from_calibration: [FAIL] Phi disagrees by rel %.2e -- too large to be ' ...
           'rounding.\n' ...
           '    Phi_implied = fp.r*(1-fp.tau_K)            = %.16g\n' ...
           '    Phi_struct  = rho + g_N + sigma*g_A        = %.16g\n' ...
           '  The stored sigma is not the Euler residual of the stored ' ...
           '(r, tau_K, rho, g_N, g_A). Run check_euler_phi.m on this file ' ...
           'for the full tiered diagnostic (including the no-g_N test) ' ...
           'before relying on it.'], rel_gap, Phi_implied, Phi_struct);
end

Phi = Phi_implied;

end
