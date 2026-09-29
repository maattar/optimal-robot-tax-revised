function rep = compute_policy_report_pol(x, cfg, cp)
% COMPUTE_POLICY_REPORT_POL  A full BGP[ and welfare report at any
% point in the instrument space, shaped by pol_config.
%
%   rep = compute_policy_report_pol(x, cfg, cp)
%
%   x   the free instrument values, in the order cfg.free lists them.
%   cp  as in eval_policy_pol. 
%









params = cp;
for i = 1:numel(cfg.free)
    params.(cfg.free{i}) = x(i);
end
params.verbose = false;
params.s_n     = 1 - params.s_a1 - params.s_a2;

% POL instruments
tau_R_req = params.tau_R;
tau_K_req = params.tau_K;
tau_L_req = params.tau_L;
w_min_req = params.w_min;
sigma_req = params.sigma;   
[sol_comp, sol_mid, regime, bgp_ok, params] = ...
    solve_policy_bgp_pol(params, cfg.needs_euler);

if ~bgp_ok
    rep = failed_report(tau_R_req, tau_K_req, tau_L_req, w_min_req, sigma_req, regime, cfg);
    return
end

sc   = compute_scalars(params);
welf = compute_welfare_regimes(sol_comp, sol_mid, regime, sc, params);

if ~welf.converged
    rep = failed_report(tau_R_req, tau_K_req, tau_L_req, w_min_req, sigma_req, regime, cfg);
    return
end







s_a1 = params.s_a1;  s_a2 = params.s_a2;  s_n = params.s_n;
u_a1 = welf.u_a1;    u_a2 = welf.u_a2;
X    = welf.X;       kx   = welf.kx;      y  = welf.y;

L_a1 = s_a1 - u_a1;
L_a2 = s_a2 - u_a2;
k    = kx - X;

if y > 0
    K2Y         = kx / y;
    labor_share = (welf.w_a1 * L_a1 + welf.w_a2 * L_a2 + s_n * welf.w_n) / y;
else
    K2Y = NaN; labor_share = NaN;
end

if kx > 0
    IRshare = X / kx;
else
    IRshare = NaN;
end
redist_budget = welf.T_AN - (u_a1 + u_a2) * welf.b;

%  automated share
switch regime
    case 'competitive'
        task_regime = sol_comp.regime;
        w_a1_eff    = sol_comp.w_a1;
        w_a2_eff    = sol_comp.w_a2;
    case 'mw_binding'
        task_regime = sol_mid.regime;
        w_a1_eff    = sol_mid.w_min;
        w_a2_eff    = sol_mid.w_a2;
    case 'joint_binding'
        task_regime = sol_mid.regime;
        w_a1_eff    = sol_mid.w_min;
        w_a2_eff    = sol_mid.w_a2;
    otherwise
        task_regime = 'failed';
        w_a1_eff    = NaN;
        w_a2_eff    = NaN;
end

robot_area = local_task_area(w_a1_eff, w_a2_eff, task_regime, sc);
total_area = sc.lam_bar * (sc.eta_bar - 1);

if total_area > 0 && isfinite(robot_area)
    automation_share = robot_area / total_area;
else
    automation_share = NaN;
end

%%avg hand-to-mouth income
n_En = 0.4 - s_a1 - s_a2;
avg_HtM = (u_a1 * welf.c_Ua1 + u_a2 * welf.c_Ua2 + ...
           L_a1 * welf.c_Ea1 + L_a2 * welf.c_Ea2 + ...
           n_En * welf.c_En) / 0.4;

%  AVERAGe Ricardian income, before saving
gro = params.g_A + params.g_N;
csh_Q3 = 2.65164070268479  / 100;
csh_Q4 = 8.71726881007623  / 100;
csh_Q5 = 88.631090487239   / 100;
kx_Q3 = (csh_Q3 / 0.2) * kx;
kx_Q4 = (csh_Q4 / 0.2) * kx;
kx_Q5 = (csh_Q5 / 0.2) * kx;
ya_Q3 = welf.c_Q3 + gro * kx_Q3;
ya_Q4 = welf.c_Q4 + gro * kx_Q4;
ya_Q5 = welf.c_Q5 + gro * kx_Q5;
avg_Ricardian = (ya_Q3 + ya_Q4 + ya_Q5) / 3;


rep = welf;
rep.tau_R            = tau_R_req;
rep.tau_K            = tau_K_req;
rep.tau_L            = tau_L_req;
rep.w_min            = w_min_req;
rep.sigma            = sigma_req;
rep.k                = k;
rep.L_a1             = L_a1;
rep.L_a2             = L_a2;
rep.K2Y              = K2Y;
rep.labor_share      = labor_share;
rep.IRshare          = IRshare;
rep.redist_budget    = redist_budget;
rep.task_regime      = task_regime;
rep.automation_share = automation_share;
rep.avg_HtM          = avg_HtM;
rep.avg_Ricardian    = avg_Ricardian;
rep.converged        = true;

if cfg.needs_euler
    r_used                = params.r;
    rep.r                 = r_used;
    rep.r_plus_delta      = r_used + params.delta;
    rep.after_tax_return  = (1 - tau_K_req) * r_used;   % equals Phi
    rep.tau_K_revenue     = tau_K_req * r_used * kx;
end

end  % main function


% The area of the robot region in the task space.
function area = local_task_area(w_a1, w_a2, task_regime, sc)
    if any(isnan([w_a1, w_a2]))
        area = NaN;
        return
    end
    p = sc.p; eta_bar = sc.eta_bar; lam_bar = sc.lam_bar;
    lam_hat = p / w_a1;
    switch task_regime
        case 'no_auto'
            area = 0;
        case 'case1'
            area = (w_a2 / p) * (lam_bar^2 - lam_hat^2) / 2 - (lam_bar - lam_hat);
        case 'case2'
            lam_prime = (p / w_a2) * eta_bar;
            area = (w_a2 / p) * (lam_prime^2 - lam_hat^2) / 2 - (lam_prime - lam_hat) ...
                 + (eta_bar - 1) * (lam_bar - lam_prime);
        otherwise
            area = NaN;
    end
end




function rep = failed_report(tau_R, tau_K, tau_L, w_min, sigma, regime, cfg)
    rep.SWF   = -Inf;
    rep.regime = regime;   
    rep.tau_R = tau_R;  rep.tau_K = tau_K;  rep.tau_L = tau_L;  rep.w_min = w_min;
    
    rep.sigma = sigma;
    rep.u_a1 = NaN; rep.u_a2 = NaN; rep.w_a1 = NaN; rep.w_a2 = NaN; rep.w_n = NaN;
    rep.d = NaN; rep.b = NaN; rep.y = NaN; rep.X = NaN; rep.kx = NaN; rep.T_AN = NaN;
    rep.c_Ua1 = NaN; rep.c_Ua2 = NaN; rep.c_Ea1 = NaN; rep.c_Ea2 = NaN; rep.c_En = NaN;
    rep.c_Q3 = NaN; rep.c_Q4 = NaN; rep.c_Q5 = NaN;
    rep.k = NaN; rep.L_a1 = NaN; rep.L_a2 = NaN; rep.K2Y = NaN; rep.labor_share = NaN;
    rep.IRshare = NaN; rep.redist_budget = NaN; rep.task_regime = 'failed';
    rep.automation_share = NaN; rep.avg_HtM = NaN; rep.avg_Ricardian = NaN;
    if cfg.needs_euler
        rep.r = NaN; rep.r_plus_delta = NaN; rep.after_tax_return = NaN;
        rep.tau_K_revenue = NaN;
    end
    rep.converged = false;
end
