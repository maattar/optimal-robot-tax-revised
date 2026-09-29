function [sol_comp, sol_mw, bgp_ok] = solve_bgp(params)
% SOLVE_BGP  Competitive BGP, binding check, then min-wage BGP.
%
%   [sol_comp, sol_mw, bgp_ok] = solve_bgp(params)
%
%
%   bgp_ok is false, and sol_mw.regime is set to 'penalty', if any of
%   three things goes wrong: the competitive BGP fails to converge or to
%   confirm its regime; the floor fails to bind, meaning w_min does not
%   strictly exceed the competitive type-1 wage; or the min-wage BGP
%   fails, which also covers the joint-binding region where it cannot
%   bracket w_a2.  The objective must test bgp_ok before touching the
%   moments and penalise the draw when it is false.
%
%   params.w_min is supplied by the caller; it is one of the free
%   calibration parameters.

params.tau_R   = 0;
params.s_n     = 1 - params.s_a1 - params.s_a2;
params.verbose = false;

w_min = params.w_min;

% --- Competitive BGP ------------------------------------------------------
sc       = compute_scalars(params);
sol_comp = solve_competitive_bgp(sc, params);

if ~sol_comp.converged || ~sol_comp.regime_ok
    sol_mw         = struct();
    sol_mw.regime  = 'penalty';
    sol_mw.converged = false;
    sol_mw.u_a1    = NaN;
    bgp_ok         = false;
    return
end

% --- The floor has to bind ------------------------------------------------
% If it does not, the competitive BGP already clears the type-1 market and
% there is no unemployment: the wrong variant of the model for this
% benchmark, not a solution to be scored.
if w_min <= sol_comp.w_a1
    sol_mw         = struct();
    sol_mw.regime  = 'penalty';
    sol_mw.converged = false;
    sol_mw.u_a1    = NaN;
    bgp_ok         = false;
    return
end

% --- Min-wage BGP ---------------------------------------------------------
sol_mw = solve_minwage_bgp(params, w_min);

if ~sol_mw.converged || strcmp(sol_mw.regime, 'failed')
    sol_mw.regime = 'penalty';
    bgp_ok        = false;
    return
end

bgp_ok = true;

end
