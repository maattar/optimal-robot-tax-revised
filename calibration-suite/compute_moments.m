function moms = compute_moments(sol_comp, sol_mw, sc, params)
% COMPUTE_MOMENTS  Model moments and welfare at a solved BGP.
%
%   moms = compute_moments(sol_comp, sol_mw, sc, params)
%
%   Returns every moment the calibration might use, targeted or not, plus
%   the utilitarian social welfare function. All quantities come from the
%   min-wage BGP; sol_comp is read only to confirm convergence.
%
%   params supplies the fiscal and preference parameters that are not
%   carried in sc: tau_K, tau_L, nu, sigma, g_A, g_N.
%
%   Every field is NaN if either solver failed, or if the solution implies
%   a negative population mass, no investment, or non-positive consumption
%   for some group.

if isempty(sol_mw) || ~sol_mw.converged || ~sol_comp.converged
    moms = nan_moms();
    return
end

y    = sol_mw.y;
k    = sol_mw.k;
X    = sol_mw.X;
w_a1 = sol_mw.w_min;   % the floor is the type-1 wage in this equilibrium
w_a2 = sol_mw.w_a2;
w_n  = sol_mw.w_n;
u_a1 = sol_mw.u_a1;
L_a1 = sol_mw.L_a1;

r    = sc.r;
delta = sc.delta;
s_n  = sc.s_n;
s_a1 = sc.s_a1;
s_a2 = sc.s_a2;

tau_K = params.tau_K;
tau_L = params.tau_L;
nu    = params.nu;
sigma = params.sigma;
g_A   = params.g_A;
g_N   = params.g_N;

kx   = k + X;                    % capital plus robots
n_En = 0.4 - s_a1 - s_a2;       % hand-to-mouth non-automatable workers

if n_En < 0
    moms = nan_moms();
    return
end

% --- Government budget ----------------------------------------------------
% Revenue from capital and labour income (tau_R is zero at the benchmark)
% funds the unemployment benefit and a lump-sum transfer to everyone in the
% bottom three quintiles who is working.  The mass of transfer recipients
% shrinks with unemployment, so d has to be solved for, not assumed.
T_AN = tau_K * r * kx  +  tau_L * (w_a1*L_a1 + w_a2*s_a2 + s_n*w_n);

b   = nu * w_a1;        % unemployment benefit, untaxed
chi_bar = 0.6;         % transfer-eligible share, bottom three quintiles
N_d = chi_bar - u_a1;  % mass receiving the transfer

if N_d <= 0
    moms = nan_moms();
    return
end
d = (T_AN - u_a1 * b) / N_d;

% --- Capital held by the Ricardian quintiles ------------------------------
% Aggregate shares converted to per-capita holdings by dividing by 0.2.
csh_Q3 = 2.65164070268479  / 100;
csh_Q4 = 8.71726881007623  / 100;
csh_Q5 = 88.631090487239   / 100;

kx_Q3 = (csh_Q3 / 0.2) * kx;
kx_Q4 = (csh_Q4 / 0.2) * kx;
kx_Q5 = (csh_Q5 / 0.2) * kx;

% --- Aggregate moments ----------------------------------------------------

K2Y = kx / y;

labor_share = (w_a1*L_a1 + w_a2*s_a2 + s_n*w_n) / y;

% Robot share of gross investment; both stocks grow at g_A + g_N.
i_K = (g_A + g_N + delta) * k;
i_R = (g_A + g_N + delta) * X;
if (i_K + i_R) <= 0
    moms = nan_moms();
    return
end
IRshare = i_R / (i_K + i_R);

wage_ratio = w_a2 / w_a1;

% Type-1 unemployment is the model's only source of joblessness, so it is
% also the aggregate rate.
u1 = u_a1;

% --- Incomes by group -----------------------------------------------------
% Seven groups: unemployed and employed type-1, type-2, hand-to-mouth
% non-automatable, and the three Ricardian quintiles.
%

yb_Ua1 = 0;
yb_Ea1 = w_a1;
yb_Ea2 = w_a2;
yb_En  = w_n;
yb_Q3  = w_n + r*kx_Q3;
yb_Q4  = w_n + r*kx_Q4;
yb_Q5  = w_n + r*kx_Q5;

% After taxes and transfers.  The benefit is untaxed, and the transfer
% reaches the bottom three quintiles only.
ya_Ua1 = b;
ya_Ea1 = (1-tau_L)*w_a1 + d;
ya_Ea2 = (1-tau_L)*w_a2 + d;
ya_En  = (1-tau_L)*w_n  + d;
ya_Q3  = (1-tau_L)*w_n  + (1-tau_K)*r*kx_Q3 + d;
ya_Q4  = (1-tau_L)*w_n  + (1-tau_K)*r*kx_Q4;
ya_Q5  = (1-tau_L)*w_n  + (1-tau_K)*r*kx_Q5;

% --- Consumption ----------------------------------------------------------
% Hand-to-mouth groups consume their income; the Ricardian quintiles save
% enough to keep their capital growing at the BGP rate.
c_Ua1 = ya_Ua1;
c_Ea1 = ya_Ea1;
c_Ea2 = ya_Ea2;
c_En  = ya_En;
c_Q3  = ya_Q3 - (g_A + g_N) * kx_Q3;
c_Q4  = ya_Q4 - (g_A + g_N) * kx_Q4;
c_Q5  = ya_Q5 - (g_A + g_N) * kx_Q5;

% CRRA utility is undefined at zero consumption.
if any([c_Ua1, c_Ea1, c_Ea2, c_En, c_Q3, c_Q4, c_Q5] <= 0)
    moms = nan_moms();
    return
end

% --- Welfare --------------------------------------------------------------
n_Ua1 = u_a1;
n_Ea1 = s_a1 - u_a1;
n_Ea2 = s_a2;

crra = @(c, n) n * (c.^(1-sigma) - 1) / (1 - sigma);

SWF = crra(c_Ua1, n_Ua1) + crra(c_Ea1, n_Ea1) + crra(c_Ea2, n_Ea2) + ...
      crra(c_En,  n_En)  + crra(c_Q3,  0.2)   + crra(c_Q4,  0.2)   + ...
      crra(c_Q5,  0.2);

% --- Distribution ---------------------------------------------------------

n_vec  = [u_a1; n_Ea1; n_Ea2; n_En; 0.2; 0.2; 0.2];
yb_vec = [yb_Ua1; yb_Ea1; yb_Ea2; yb_En; yb_Q3; yb_Q4; yb_Q5];
ya_vec = [ya_Ua1; ya_Ea1; ya_Ea2; ya_En; ya_Q3; ya_Q4; ya_Q5];

gini_before = weighted_gini(yb_vec, n_vec);
gini_after = weighted_gini(ya_vec, n_vec);

bq_share_after = bottom_income_share(ya_vec, n_vec, 0.20);
bq_share_before = bottom_income_share(yb_vec, n_vec, 0.20);

top20_share_before = top_income_share(yb_vec, n_vec, 0.20);
top20_share_after  = top_income_share(ya_vec, n_vec, 0.20);

% --- Wages ----------------------------------------------------------------
% Three wage levels among the employed: the floor with mass L_a1, the
% type-2 wage with mass s_a2, and w_n for every non-automatable worker,
% hand-to-mouth and Ricardian alike.
w_vals  = [w_a1; w_a2; w_n];
w_wts   = [L_a1; s_a2; s_n];
p50_p10 = wage_percentile_ratio(w_vals, w_wts, 50, 10);

total_emp      = 1 - u_a1;
mean_wage      = (w_a1*L_a1 + w_a2*s_a2 + w_n*s_n) / total_emp;
min_mean_ratio = w_a1 / mean_wage;

emp_minwage_sh = L_a1;

% --- Pack -----------------------------------------------------------------
moms.K2Y          = K2Y;
moms.labor_share  = labor_share;
moms.IRshare      = IRshare;
moms.wage_ratio   = wage_ratio;
moms.u_a1         = u1;

moms.min_mean_ratio  = min_mean_ratio;
moms.gini_before     = gini_before;
moms.gini_after      = gini_after;
moms.bq_share_after  = bq_share_after;
moms.bq_share_before = bq_share_before;
moms.top20_share_before = top20_share_before;
moms.top20_share_after  = top20_share_after;
moms.p50_p10         = p50_p10;
moms.emp_minwage_sh  = emp_minwage_sh;

moms.SWF = SWF;

end


function moms = nan_moms()
    moms.K2Y            = NaN;
    moms.labor_share    = NaN;
    moms.IRshare        = NaN;
    moms.wage_ratio     = NaN;
    moms.u_a1           = NaN;
    moms.min_mean_ratio = NaN;
    moms.gini_before    = NaN;
    moms.gini_after     = NaN;
    moms.bq_share_after = NaN;
    moms.bq_share_before    = NaN;
    moms.top20_share_before = NaN;
    moms.top20_share_after  = NaN;
    moms.p50_p10        = NaN;
    moms.emp_minwage_sh = NaN;
    moms.SWF            = NaN;
end


function G = weighted_gini(y, n)
% Gini for a discrete grouped distribution.

    [y_s, idx] = sort(y);
    n_s  = n(idx);
    n_tot = sum(n_s);
    y_tot = sum(y_s .* n_s);

    if y_tot <= 0
        G = 0;
        return
    end

    F = [0; cumsum(n_s) / n_tot];         % cumulative population share
    L = [0; cumsum(y_s .* n_s) / y_tot];  % cumulative income share

    G = 1 - sum( (F(2:end) - F(1:end-1)) .* (L(2:end) + L(1:end-1)) );
end


function share = bottom_income_share(y, n, frac)

    [y_s, idx] = sort(y);
    n_s  = n(idx);
    n_tot = sum(n_s);
    y_tot = sum(y_s .* n_s);

    if y_tot <= 0
        share = 0;
        return
    end

    target     = frac * n_tot;
    income_bot = 0;
    pop_bot    = 0;

    for i = 1:numel(y_s)
        remaining = target - pop_bot;
        if n_s(i) <= remaining
            income_bot = income_bot + y_s(i) * n_s(i);
            pop_bot    = pop_bot    + n_s(i);
        else
            income_bot = income_bot + y_s(i) * remaining;
            break
        end
        if pop_bot >= target
            break
        end
    end

    share = income_bot / y_tot;
end


function share = top_income_share(y, n, frac)
    share = 1 - bottom_income_share(y, n, 1 - frac);
end


function ratio = wage_percentile_ratio(w_vals, w_wts, p_high, p_low)

    [w_s, idx] = sort(w_vals);
    wt_s  = w_wts(idx);
    wt_tot = sum(wt_s);
    cum_wt = cumsum(wt_s) / wt_tot;

    w_high = find_pctile(w_s, cum_wt, p_high / 100);
    w_low  = find_pctile(w_s, cum_wt, p_low  / 100);

    if w_low <= 0
        ratio = NaN;
    else
        ratio = w_high / w_low;
    end
end

function wp = find_pctile(w_s, cum_wt, frac)
    idx = find(cum_wt >= frac, 1, 'first');
    if isempty(idx)
        wp = w_s(end);
    else
        wp = w_s(idx);
    end
end
