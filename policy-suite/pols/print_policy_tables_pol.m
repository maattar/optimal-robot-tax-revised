function [T1, T2, cev] = print_policy_tables_pol(rep0, rep1, cfg, cp, pol_id, title_suffix)
% PRINT_POLICY_TABLES_POL  The two comparison tables and the aggregate CEV
% block, in the same form for every policy.
%
%   [T1, T2, cev] = print_policy_tables_pol(rep0, rep1, cfg, cp, pol_id, title_suffix)
%
% Table 1 sets the status quo against the optimum in aggregate; Table 2
% gives consumption levels and CEV by worker group; the aggregate CEV
% follows. 
if ~rep0.converged
    error(['print_policy_tables_pol: status-quo report has not converged ' ...
           '(regime=%s) -- nothing to print. Caller should have checked ' ...
           'this before calling.'], rep0.regime);
end
if ~rep1.converged
    error(['print_policy_tables_pol: optimal report has not converged ' ...
           '(regime=%s) -- nothing to print. Caller should have checked ' ...
           'this before calling.'], rep1.regime);
end
if nargin < 6, title_suffix = ''; end


Group = {'Unemployed, type-1, Hand-to-Mouth (U_{a1})'; ...
         'Unemployed, type-2, Hand-to-Mouth (U_{a2})'; ...
         'Employed, type-1, Hand-to-Mouth (E_{a1})'; ...
         'Employed, type-2, Hand-to-Mouth (E_{a2})'; ...
         'Employed, type-n, Hand-to-Mouth (E_n)'; ...
         'Ricardian, Q3'; 'Ricardian, Q4'; 'Ricardian, Q5'};

c_bench = [rep0.c_Ua1; rep0.c_Ua2; rep0.c_Ea1; rep0.c_Ea2; rep0.c_En; rep0.c_Q3; rep0.c_Q4; rep0.c_Q5];
c_opt   = [rep1.c_Ua1; rep1.c_Ua2; rep1.c_Ea1; rep1.c_Ea2; rep1.c_En; rep1.c_Q3; rep1.c_Q4; rep1.c_Q5];
PopShare_StatusQuo = [rep0.u_a1; rep0.u_a2; rep0.L_a1; rep0.L_a2; 0.4-cp.s_a1-cp.s_a2; 0.2; 0.2; 0.2];
PopShare_Optimal   = [rep1.u_a1; rep1.u_a2; rep1.L_a1; rep1.L_a2; 0.4-cp.s_a1-cp.s_a2; 0.2; 0.2; 0.2];

cev = compute_cev_pol(rep0, rep1, PopShare_StatusQuo, PopShare_Optimal);

% =========================================================================
% Table 1: the status quo against the optimum, in aggregate
% =========================================================================
tauR_label = 'Robot Tax  \tau_R';
if ~any(strcmp('tau_R', cfg.free)), tauR_label = [tauR_label '  (FIXED)']; end
Variable   = {tauR_label};
level_vars = {rep0.tau_R}; opt_vars = {rep1.tau_R};

wmin_label = 'Wage Floor  w_{min}';
if ~any(strcmp('w_min', cfg.free)), wmin_label = [wmin_label '  (FIXED)']; end
if cfg.needs_euler
    Variable = [Variable; {'Capital Tax  \tau_K'; 'Labor Tax  \tau_L'; wmin_label}];
    level_vars = [level_vars, {rep0.tau_K, rep0.tau_L, rep0.w_min}];
    opt_vars   = [opt_vars,   {rep1.tau_K, rep1.tau_L, rep1.w_min}];
else
    Variable = [Variable; {wmin_label}];
    level_vars = [level_vars, {rep0.w_min}];
    opt_vars   = [opt_vars,   {rep1.w_min}];
end
n_level = numel(Variable);

Variable   = [Variable; {'Type-1 Wage  w_{a1}'; 'Type-2 Wage  w_{a2}'}];
level_vars = [level_vars, {rep0.w_a1, rep0.w_a2}];
opt_vars   = [opt_vars,   {rep1.w_a1, rep1.w_a2}];
n_level    = n_level + 2;

if cfg.needs_euler
    Variable   = [Variable; {'Interest Rate  r  (ENDOGENOUS)'; 'After-Tax Return  (1-\tau_K)r'}];
    level_vars = [level_vars, {rep0.r, rep0.after_tax_return}];
    opt_vars   = [opt_vars,   {rep1.r, rep1.after_tax_return}];
    n_level    = n_level + 2;
end

Variable = [Variable; {'Social Welfare'; 'Fraction of Automated Tasks'; ...
            'Type-1 Unemployment Rate'; 'Type-2 Unemployment Rate'; ...
            'Labor Income Share'; 'Capital-Output Ratio'; ...
            'Output (per effective worker)'; 'Avg. Hand-to-Mouth Income'; ...
            'Avg. Ricardian Income'; 'Redistribution Budget'}];
rest_bench = [rep0.SWF, rep0.automation_share, rep0.u_a1, rep0.u_a2, ...
              rep0.labor_share, rep0.K2Y, rep0.y, rep0.avg_HtM, rep0.avg_Ricardian, rep0.redist_budget];
rest_opt   = [rep1.SWF, rep1.automation_share, rep1.u_a1, rep1.u_a2, ...
              rep1.labor_share, rep1.K2Y, rep1.y, rep1.avg_HtM, rep1.avg_Ricardian, rep1.redist_budget];

Benchmark = [cell2mat(level_vars)'; rest_bench'];
Optimal   = [cell2mat(opt_vars)';   rest_opt'];

is_level = false(numel(Variable), 1); is_level(1:n_level) = true;
is_pp    = false(numel(Variable), 1);
is_pp(n_level + 3) = true; is_pp(n_level + 4) = true;   % u_a1, u_a2 rows

Change = cell(numel(Variable), 1);
for i = 1:numel(Variable)
    if is_level(i)
        Change{i} = sprintf('%+.6f  (level)', Optimal(i) - Benchmark(i));
    else
        Change{i} = format_change(Benchmark(i), Optimal(i), is_pp(i));
    end
end

Variable  = [Variable;  {'Aggregate CEV (%)'}];
Benchmark = [Benchmark; 0];
Optimal   = [Optimal;   cev.CEV_agg_pct];
Change    = [Change;    {'reallocation gain -- see the Aggregate CEV block'}];

T1 = make_table_or_struct(Variable, Benchmark, Optimal, Change, ...
    {'Variable', 'Benchmark', 'Optimal', 'Change'});

fprintf('\n===========================================================\n');
if isempty(title_suffix)
    fprintf('  TABLE 1: STATUS QUO vs. OPTIMAL %s POLICY  [%s]\n', ...
            upper(pol_id), cfg.manuscript_label);
else
    fprintf('  TABLE 1: STATUS QUO vs. OPTIMAL %s POLICY  [%s]  (%s)\n', ...
            upper(pol_id), cfg.manuscript_label, title_suffix);
end
fprintf('===========================================================\n');

fprintf('  Regime:  bench = %-14s   opt = %s\n', rep0.regime, rep1.regime);

print_table(Variable, Benchmark, Optimal, Change);

BLOCKED_EPS = 1e-9;
fprintf('  %-38s  %12s  %12s\n', 'Automation Blocked', ...
        blocked_str(rep0.automation_share, BLOCKED_EPS), ...
        blocked_str(rep1.automation_share, BLOCKED_EPS));
if cfg.needs_euler
    fprintf(['  Reading note: the after-tax return is invariant at Phi by\n' ...
             '  construction. A zero change there is CORRECT, not a stuck\n' ...
             '  variable.\n']);
end

CEV_pct = cev.CEV_g_pct;

T2 = make_table_or_struct(Group, PopShare_StatusQuo, PopShare_Optimal, ...
    c_bench, c_opt, CEV_pct, ...
    {'Group', 'PopShare_StatusQuo', 'PopShare_Optimal', ...
     'Consumption_StatusQuo', 'Consumption_Optimal', 'CEV_pct'});

fprintf('\n===========================================================\n');
fprintf('  TABLE 2: CONSUMPTION LEVELS AND CEV BY GROUP\n');
fprintf('===========================================================\n');

fprintf('  %-46s  %10s  %10s  %10s  %10s  %10s\n', ...
        'Group', 'pop0', 'pop1', 'c0', 'c1', 'CEV');
for i = 1:numel(Group)
    if isnan(CEV_pct(i))
        cev_str = '       N/A';
    else
        cev_str = sprintf('%+9.3f%%', CEV_pct(i));
    end
    fprintf('  %-46s  %10.4f  %10.4f  %10.6f  %10.6f  %10s\n', ...
            Group{i}, PopShare_StatusQuo(i), PopShare_Optimal(i), ...
            c_bench(i), c_opt(i), cev_str);
end
fprintf('\n');

MASS_TOL     = 1e-8;
IDENTITY_TOL = 1e-9;

fprintf('===========================================================\n');
fprintf('  AGGREGATE CEV\n');
fprintf('===========================================================\n');
if ~cev.converged
    fprintf(['  CEV^agg:  n/a -- the consumption vectors are not usable for a\n' ...
             '  welfare comparison (a report failed to converge, or a populated\n' ...
             '  group has non-positive consumption).\n\n']);
else
    fprintf('  %-30s  %+14.6f%%\n', 'CEV^agg', cev.CEV_agg_pct);
    fprintf('  %-30s  %14.6f\n\n', 'sigma', cev.sigma);
    fprintf('  %-30s  %16s  %16s\n', 'Component', 'bench', 'opt');
    fprintf('  %-30s  %16.10f  %16.10f\n', 'sum_g m_g c_g^(1-sigma)', cev.S_sq, cev.S_opt);
    fprintf('  %-30s  %16.10f  %16.10f\n', 'active population mass', cev.mass_sq, cev.mass_opt);
    fprintf('  %-30s  %16.10f  %16.10f\n', 'SWF', rep0.SWF, rep1.SWF);

    if abs(cev.mass_sq - 1) > MASS_TOL || abs(cev.mass_opt - 1) > MASS_TOL
        fprintf(['\n  WARNING: the active population shares do not sum to 1 ' ...
                 '(bench %.12f, opt %.12f).\n  CEV^agg assumes they do -- ' ...
                 'check the group weights before quoting this number.\n'], ...
                cev.mass_sq, cev.mass_opt);
    end

    if isfinite(cev.CEV_agg_from_SWF)
        fprintf('\n  Consistency check (exact identity, S = 1 + (1-sigma)*SWF):\n');
        fprintf('  %-30s  %+14.6f%%\n', 'CEV^agg recovered from SWF', 100*cev.CEV_agg_from_SWF);
        if cev.identity_gap > IDENTITY_TOL
            fprintf(['  WARNING: gap = %.3e, larger than rounding. The consumption\n' ...
                     '  levels and the SWF disagree -- do not quote CEV^agg until\n' ...
                     '  that is resolved.\n'], cev.identity_gap);
        else
            fprintf('  %-30s  %14.3e  (OK)\n', 'gap', cev.identity_gap);
        end
    end

    fprintf(['\n  Reading note: CEV^agg credits workers MOVING between groups ' ...
             '(e.g.\n  unemployed into employment) on top of the change in what each\n' ...
             '  group consumes. The CEV_g column above holds the composition of the\n' ...
             '  population fixed and so answers a different question; the two need\n' ...
             '  not agree in sign, and CEV^agg is not an average of them.\n\n']);
end

end  


function T = make_table_or_struct(varargin)
    names = varargin{end};
    cols  = varargin(1:end-1);
    try
        T = table(cols{:}, 'VariableNames', names);
        return
    catch
        
    end
    T = struct();
    for i = 1:numel(names)
        T.(names{i}) = cols{i};
    end
end

function print_table(Variable, Benchmark, Optimal, Change)
  
    fprintf('  %-38s  %12s  %12s  %s\n', 'Variable', 'bench', 'opt', 'Change');
    for i = 1:numel(Variable)
        fprintf('  %-38s  %12.6f  %12.6f  %s\n', ...
                Variable{i}, Benchmark(i), Optimal(i), Change{i});
    end
end

function s = format_change(v0, v1, is_pp)
    if is_pp
        d = 100 * (v1 - v0);
        s = sprintf('%.2f p.p. %s', abs(d), word_for(d));
    else
        if v0 == 0 || ~isfinite(v0)
            s = 'n/a (benchmark = 0)';
            return
        end
        d = 100 * (v1 - v0) / abs(v0);
        s = sprintf('%.2f%% %s', abs(d), word_for(d));
    end
end

function w = word_for(d)
    if d >= 0, w = 'increase'; else, w = 'decrease'; end
end

function s = blocked_str(automation_share, eps_tol)
    if automation_share <= eps_tol
        s = 'Yes';
    else
        s = 'No';
    end
end
