function [sol_comp, sol_mid, regime, bgp_ok] = solve_policy_bgp_regimes(params)
% SOLVE_POLICY_BGP_REGIMES
%   regime   competitive, mw_binding, joint_binding or failed.
%   bgp_ok   true for any of the three that is not failure.
params.s_n     = 1 - params.s_a1 - params.s_a2;
params.verbose = false;

w_min = params.w_min;

% competitive path solved first
sc       = compute_scalars(params);
sol_comp = solve_competitive_bgp(sc, params);

if ~sol_comp.converged || ~sol_comp.regime_ok
    sol_mid = make_failed_mid();
    regime  = 'failed';
    bgp_ok  = false;
    return
end

sol_comp = recover_post_solution(sol_comp, sc, params);   % adds output, capital, the wage and robot input

if ~sol_comp.checks.all_ok
    sol_mid = make_failed_mid();
    regime  = 'failed';
    bgp_ok  = false;
    return
end

w_a1_c = sol_comp.w_a1;

% regime check
if w_min <= w_a1_c
    sol_mid = struct('regime', 'competitive', 'converged', true, ...
                      'u_a1', 0.0, 'u_a2', 0.0);
    regime  = 'competitive';
    bgp_ok  = true;
    return
end



sol_try2 = solve_minwage_bgp(params, w_min);

if sol_try2.converged && ~strcmp(sol_try2.regime, 'failed')
    sol_try2.u_a2 = 0.0;              % tHIS market clears
    sol_try2.L_a2 = params.s_a2;
    sol_mid = sol_try2;
    regime  = 'mw_binding';
    bgp_ok  = true;
    return
end



sol_try3 = solve_joint_binding_bgp(params, w_min);

if sol_try3.converged && ~strcmp(sol_try3.regime, 'failed')
    sol_mid = sol_try3;
    regime  = 'joint_binding';
    bgp_ok  = true;
    return
end

% FFailure.
sol_mid = make_failed_mid();
regime  = 'failed';
bgp_ok  = false;

end  



function s = make_failed_mid()
    s.regime    = 'failed';
    s.converged = false;
    s.u_a1      = NaN;
    s.u_a2      = NaN;
end
