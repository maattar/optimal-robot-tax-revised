function welf = eval_policy_pol(x, cfg, cp)
% EVAL_POLICY_POL  Evaluates one point in the instrument space. The same
% function serves every policy, taking its shape from pol_config.



if numel(x) ~= numel(cfg.free)
    error(['eval_policy_pol: x has %d element(s) but cfg (%s) expects ' ...
           '%d free instrument(s): %s.'], numel(x), cfg.id, ...
          numel(cfg.free), strjoin(cfg.free, ', '));
end

params = cp;
for i = 1:numel(cfg.free)
    params.(cfg.free{i}) = x(i);
end

[sol_comp, sol_mid, regime, bgp_ok, params] = ...
    solve_policy_bgp_pol(params, cfg.needs_euler);

if ~bgp_ok
    welf = struct( ...
        'converged', false,   'SWF', -Inf,  'regime', regime, ...
        'tau_R', params.tau_R, 'tau_K', params.tau_K, ...
        'tau_L', params.tau_L, 'w_min', params.w_min, ...
        'r_used', NaN, ...
        'u_a1', NaN,          'u_a2', NaN,   ...
        'w_a1', NaN,          'w_a2',  NaN,  'w_n',  NaN,   ...
        'd',    NaN,          'b',     NaN,  'y',    NaN,   ...
        'X',    NaN,          'kx',    NaN,  'T_AN', NaN,   ...
        'c_Ua1', NaN, 'c_Ua2', NaN, 'c_Ea1', NaN, 'c_Ea2', NaN, ...
        'c_En', NaN, 'c_Q3',  NaN, 'c_Q4',  NaN, 'c_Q5',  NaN);
    return
end


sc   = compute_scalars(params);
welf = compute_welfare_regimes(sol_comp, sol_mid, regime, sc, params);

welf.tau_K  = params.tau_K;
welf.tau_L  = params.tau_L;
welf.r_used = params.r;

end
