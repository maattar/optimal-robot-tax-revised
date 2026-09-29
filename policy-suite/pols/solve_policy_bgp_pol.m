function [sol_comp, sol_mid, regime, bgp_ok, params_out] = solve_policy_bgp_pol(params, needs_euler)
% SOLVE_POLICY_BGP_POL  Settles the interest rate, then hands the balanced
% path to the regime solver. Whether the Euler equation fires is read from
% pol_config, so one function covers the policies that hold r fixed and
% those that derive it.
%
%   [sol_comp, sol_mid, regime, bgp_ok, params_out] = ...
%       solve_policy_bgp_pol(params, needs_euler)
%
% With needs_euler false (pol1, pol6) params.r is required and passes
% through untouched: tau_K never moves under those policies, so the Euler
% equation would return the calibrated rate anyway.
%
% With needs_euler true (pol2, pol3, pol4) params.Phi and
% params.tau_K are required, and any r already in params is ignored and
% overwritten with euler_closure(Phi, tau_K). Should that come back
% non-finite, tau_K having reached the pole, the regime returns as singular
% rather than raising an error -- singular is not failed: nothing has gone
% wrong, the point simply is not an economy, and a parallel sweep can step
% over it.
%
% This file adds nothing to the three-regime dispatch itself, which belongs
% to solve_policy_bgp_regimes; it only decides the interest rate first.
%
% params_out always carries the rate actually used, passed through or
% derived. Give compute_scalars and compute_welfare_regimes that struct
% rather than your own copy, which would otherwise disagree about r
% wherever the Euler step fired.
%
%   params  the structural and fiscal fields, plus Phi where the Euler
%           step is needed or r where it is not.
if nargin < 2
    error('solve_policy_bgp_pol: needs_euler is required (no silent default).');
end

params.s_n     = 1 - params.s_a1 - params.s_a2;
params.verbose = false;

if needs_euler
    need = {'Phi', 'tau_K'};
    missing = need(~cellfun(@(fn) isfield(params, fn), need));
    if ~isempty(missing)
        error(['solve_policy_bgp_pol: needs_euler = true but params is ' ...
               'missing required field(s): %s.'], strjoin(missing, ', '));
    end

    params.r = euler_closure(params.Phi, params.tau_K);

    if ~isfinite(params.r)
        params_out = params;
        sol_comp = struct('converged', false, 'regime_ok', false, ...
                           'regime', 'failed', 'w_a1', NaN, 'w_a2', NaN);
        sol_mid  = struct('regime', 'failed', 'converged', false, ...
                           'u_a1', NaN, 'u_a2', NaN);
        regime   = 'singular';
        bgp_ok   = false;
        return
    end
else
    if ~isfield(params, 'r')
        error(['solve_policy_bgp_pol: needs_euler = false but params has ' ...
               'no .r field. This pol does not free tau_K, so r must be ' ...
               'supplied by the caller (typically cp.r = fp.r from the ' ...
               'loaded calibration) -- it is not derived here.']);
    end
end

params_out = params;

[sol_comp, sol_mid, regime, bgp_ok] = solve_policy_bgp_regimes(params);

end
