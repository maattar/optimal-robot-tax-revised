function welf = compute_welfare_regimes(sol_comp, sol_mid, regime, sc, params)
% COMPUTE_WELFARE_REGIMES  Social welfare on the BGP, across all
% three regimes. Every policy reaches it through eval_policy_pol.
%
%   welf = compute_welfare_regimes(sol_comp, sol_mid, regime, sc, params)
%
% The eight groups, the first five hand to mouth:
%
%   Ua1, Ua2   unemployed type-1 and type-2, on the insurance benefit
%   Ea1, Ea2   employed type-1 and type-2, after-tax wage plus transfer
%   En         employed non-automatable, likewise
%   Q3         Ricardian: wage, after-tax capital income, plus transfer
%   Q4, Q5     Ricardian: wage and after-tax capital income
%
% Ua2 is populated only in joint_binding. 

if strcmp(regime, 'failed') || ~sol_comp.converged
    welf = failed_welf(regime, params);
    return
end

r      = sc.r;
delta  = sc.delta;
s_n    = sc.s_n;
s_a1   = sc.s_a1;
s_a2   = sc.s_a2;

tau_K  = params.tau_K;
tau_L  = params.tau_L;
tau_R  = params.tau_R;
nu     = params.nu;
sigma  = params.sigma;
g_A    = params.g_A;
g_N    = params.g_N;
w_min  = params.w_min;

% The Ricardian capital shares, from the data and unaffected by policy
csh_Q3 = 2.65164070268479  / 100;
csh_Q4 = 8.71726881007623  / 100;
csh_Q5 = 88.631090487239   / 100;


switch regime

    case 'competitive'
        y    = sol_comp.y;
        k    = sol_comp.k;
        X    = sol_comp.X;
        w_a1 = sol_comp.w_a1;       % the competitive wage, not the floor
        w_a2 = sol_comp.w_a2;
        w_n  = sol_comp.w_n;
        u_a1 = 0;
        u_a2 = 0;
        L_a1 = s_a1;
        L_a2 = s_a2;

    case 'mw_binding'
        if ~sol_mid.converged
            welf = failed_welf(regime, params);
            return
        end
        y    = sol_mid.y;
        k    = sol_mid.k;
        X    = sol_mid.X;
        w_a1 = sol_mid.w_min;       % at the floor
        w_a2 = sol_mid.w_a2;
        w_n  = sol_mid.w_n;
        u_a1 = sol_mid.u_a1;
        u_a2 = 0;                    % that market clears in this regime
        L_a1 = sol_mid.L_a1;
        L_a2 = s_a2;

    case 'joint_binding'
        if ~sol_mid.converged
            welf = failed_welf(regime, params);
            return
        end
        y    = sol_mid.y;
        k    = sol_mid.k;
        X    = sol_mid.X;
        w_a1 = sol_mid.w_min;       % both wages sit at the floor
        w_a2 = sol_mid.w_min;
        w_n  = sol_mid.w_n;
        u_a1 = sol_mid.u_a1;        % the whole of type-1 supply
        u_a2 = sol_mid.u_a2;
        L_a1 = sol_mid.L_a1;        % none of them employed
        L_a2 = sol_mid.L_a2;

    otherwise
        welf = failed_welf(regime, params);
        return
end

% Derived aggregates
kx   = k + X;
n_En = 0.4 - s_a1 - s_a2;          % HtM non-automatable pop. share

if n_En < 0
    welf = failed_welf(regime, params);
    return
end

% What each Ricardian quintile holds
kx_Q3 = (csh_Q3 / 0.2) * kx;
kx_Q4 = (csh_Q4 / 0.2) * kx;
kx_Q5 = (csh_Q5 / 0.2) * kx;

% Government revenue per effective worker.
T_AN = tau_K * r * kx ...
     + tau_L * (w_a1 * L_a1 + w_a2 * L_a2 + s_n * w_n) ...
     + tau_R * (r + delta) * X;

% The unemployment benefit and the transfer.
if (u_a1 + u_a2) > 0
    b = nu * w_min;
else
    b = 0;
end
chi_bar = 0.6;                     % transfer-eligible share, bottom three quintiles
N_d = chi_bar - u_a1 - u_a2;

if N_d <= 0
    welf = failed_welf(regime, params);
    return
end

d = (T_AN - (u_a1 + u_a2) * b) / N_d;

if d < 0
    welf = failed_welf(regime, params);
    return
end

% After-tax incomes, for the eight groups
ya_Ua1 = b;                                           
ya_Ua2 = b;                                           
ya_Ea1 = (1 - tau_L) * w_a1 + d;
ya_Ea2 = (1 - tau_L) * w_a2 + d;
ya_En  = (1 - tau_L) * w_n  + d;
ya_Q3  = (1 - tau_L) * w_n  + (1 - tau_K) * r * kx_Q3 + d;
ya_Q4  = (1 - tau_L) * w_n  + (1 - tau_K) * r * kx_Q4;
ya_Q5  = (1 - tau_L) * w_n  + (1 - tau_K) * r * kx_Q5;


gro = g_A + g_N;

c_Ua1 = ya_Ua1;
c_Ua2 = ya_Ua2;
c_Ea1 = ya_Ea1;
c_Ea2 = ya_Ea2;
c_En  = ya_En;
c_Q3  = ya_Q3 - gro * kx_Q3;
c_Q4  = ya_Q4 - gro * kx_Q4;
c_Q5  = ya_Q5 - gro * kx_Q5;

% Check that consumption is positive
c_vec = [c_Ua1; c_Ua2; c_Ea1; c_Ea2; c_En; c_Q3; c_Q4; c_Q5];
n_vec = [u_a1; u_a2; s_a1 - u_a1; s_a2 - u_a2; n_En; 0.2; 0.2; 0.2];
active = n_vec > 0;

if any(c_vec(active) <= 0)
    welf = failed_welf(regime, params);
    return
end

% Social welfare: the population-weighted sum of utilities over the eight
% groups.

crra = @(c, n) n .* (c .^ (1 - sigma) - 1) / (1 - sigma);

SWF = sum(crra(c_vec(active), n_vec(active)));

% output
welf.SWF       = SWF;
welf.regime    = regime;
welf.tau_R     = tau_R;
welf.w_min     = w_min;
welf.u_a1      = u_a1;
welf.u_a2      = u_a2;
welf.w_a1      = w_a1;
welf.w_a2      = w_a2;
welf.w_n       = w_n;
welf.d         = d;
welf.b         = b;
welf.y         = y;
welf.X         = X;
welf.kx        = kx;
welf.T_AN      = T_AN;
welf.c_Ua1     = c_Ua1;
welf.c_Ua2     = c_Ua2;
welf.c_Ea1     = c_Ea1;
welf.c_Ea2     = c_Ea2;
welf.c_En      = c_En;
welf.c_Q3      = c_Q3;
welf.c_Q4      = c_Q4;
welf.c_Q5      = c_Q5;
welf.converged = true;

end  % main function


function welf = failed_welf(regime, params)
    welf.SWF       = -Inf;
    welf.regime    = regime;
    welf.tau_R     = params.tau_R;
    welf.w_min     = params.w_min;
    welf.u_a1      = NaN;
    welf.u_a2      = NaN;
    welf.w_a1      = NaN;
    welf.w_a2      = NaN;
    welf.w_n       = NaN;
    welf.d         = NaN;
    welf.b         = NaN;
    welf.y         = NaN;
    welf.X         = NaN;
    welf.kx        = NaN;
    welf.T_AN      = NaN;
    welf.c_Ua1     = NaN;
    welf.c_Ua2     = NaN;
    welf.c_Ea1     = NaN;
    welf.c_Ea2     = NaN;
    welf.c_En      = NaN;
    welf.c_Q3      = NaN;
    welf.c_Q4      = NaN;
    welf.c_Q5      = NaN;
    welf.converged = false;
end
