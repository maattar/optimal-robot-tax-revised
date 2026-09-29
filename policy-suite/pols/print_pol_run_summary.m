function print_pol_run_summary(pol_id, SPEC, rep0, rep1, wmin_result)
% PRINT_POL_RUN_SUMMARY  
%
%   
if nargin < 5, wmin_result = []; end
cfg = pol_config(pol_id);
POP_EPS = 1e-9;

header = sprintf('FINAL RUN SUMMARY -- %s [%s]  (SPEC=%s)', ...
                  upper(pol_id), cfg.manuscript_label, SPEC);
bar = repmat('#', 1, max(63, numel(header) + 8));

fprintf('\n%s\n', bar);
fprintf('#  %s%s#\n', header, repmat(' ', 1, numel(bar) - numel(header) - 5));
fprintf('%s\n\n', bar);


fprintf('--- Regime --------------------------------------------------\n');
fprintf('  bench = %-16s   opt = %s\n\n', rep0.regime, rep1.regime);


fprintf('--- Instruments -----------------------------------------------\n');
fprintf('  %-28s  %14s  %14s\n', 'Instrument', 'bench', 'opt');
print_instr_row('tau_R', rep0.tau_R, rep1.tau_R, cfg);
print_instr_row('tau_K', rep0.tau_K, rep1.tau_K, cfg);
print_instr_row('tau_L', rep0.tau_L, rep1.tau_L, cfg);
print_instr_row('w_min', rep0.w_min, rep1.w_min, cfg);
fprintf('\n');


BLOCKED_EPS = 1e-9;
fprintf('--- Headline welfare and wages ------------------------------\n');
fprintf('  %-28s  %14s  %14s\n', 'Metric', 'bench', 'opt');
fprintf('  %-28s  %14.8f  %14.8f\n', 'SWF', rep0.SWF, rep1.SWF);
fprintf('  %-28s  %13.4f%%  %13.4f%%\n', 'Automated Task Share', ...
        100*rep0.automation_share, 100*rep1.automation_share);
fprintf('  %-28s  %14s  %14s\n', 'Automation Blocked', ...
        blocked_str(rep0.automation_share, BLOCKED_EPS), blocked_str(rep1.automation_share, BLOCKED_EPS));
if cfg.needs_euler && isfield(rep0, 'r') && isfield(rep1, 'r')
    fprintf('  %-28s  %14.6f  %14.6f\n', 'r (endogenous)', rep0.r, rep1.r);
end
fprintf('  %-28s  %14.6f  %14.6f\n', 'w_a1', rep0.w_a1, rep1.w_a1);
fprintf('  %-28s  %14.6f  %14.6f\n\n', 'w_a2', rep0.w_a2, rep1.w_a2);


fprintf('--- Consumption levels and CEV by group ----------------------\n');
Group = {'U_a1', 'U_a2', 'E_a1', 'E_a2', 'E_n', 'Q3', 'Q4', 'Q5'};

n_En0 = 0.4 - (rep0.u_a1 + rep0.L_a1) - (rep0.u_a2 + rep0.L_a2);
n_En1 = 0.4 - (rep1.u_a1 + rep1.L_a1) - (rep1.u_a2 + rep1.L_a2);
pop0 = [rep0.u_a1, rep0.u_a2, rep0.L_a1, rep0.L_a2, n_En0, 0.2, 0.2, 0.2];
pop1 = [rep1.u_a1, rep1.u_a2, rep1.L_a1, rep1.L_a2, n_En1, 0.2, 0.2, 0.2];
c0 = [rep0.c_Ua1, rep0.c_Ua2, rep0.c_Ea1, rep0.c_Ea2, rep0.c_En, rep0.c_Q3, rep0.c_Q4, rep0.c_Q5];
c1 = [rep1.c_Ua1, rep1.c_Ua2, rep1.c_Ea1, rep1.c_Ea2, rep1.c_En, rep1.c_Q3, rep1.c_Q4, rep1.c_Q5];


cev = compute_cev_pol(rep0, rep1);
CEV_pct = cev.CEV_g_pct;

fprintf('  %-6s  %10s  %10s  %10s  %10s  %10s\n', 'Group', 'pop0', 'pop1', 'c0', 'c1', 'CEV');
for i = 1:numel(Group)
    if ~isfinite(CEV_pct(i))
        cev_str = '       N/A';
    else
        cev_str = sprintf('%+9.3f%%', CEV_pct(i));
    end
    if isfinite(pop0(i))
        pop_str = sprintf('%10.4f  %10.4f', pop0(i), pop1(i));
    else
        pop_str = sprintf('%10s  %10s', 'n/a', 'n/a');
    end
    fprintf('  %-6s  %s  %10.6f  %10.6f  %10s\n', Group{i}, pop_str, c0(i), c1(i), cev_str);
end
fprintf('\n');


fprintf('--- Aggregate CEV -------------------------------------------\n');
if cev.converged
    fprintf('  %-28s  %+14.6f%%\n', 'CEV^agg', cev.CEV_agg_pct);
    fprintf('  %-28s  %14.6f\n', 'sigma', cev.sigma);
    if isfinite(cev.identity_gap) && cev.identity_gap > 1e-9
        fprintf(['  WARNING: the SWF-identity cross-check is off by %.3e. See the\n' ...
                 '  Aggregate CEV block in Table 2''s output for the components.\n'], ...
                cev.identity_gap);
    end
else
    fprintf('  %-28s  %15s\n', 'CEV^agg', 'N/A');
end
fprintf('\n');



if ~isempty(wmin_result)
    fprintf('--- w_min identification diagnostic ---------------------------\n');
    fprintf('  raw search value = %.6f   threshold (competitive w_a1) = %.6f\n', ...
            wmin_result.w_min_star, wmin_result.threshold_wage);
    if wmin_result.is_indeterminate && abs(rep1.w_min - wmin_result.w_min_star) > POP_EPS
        fprintf('  Tables above report the recommended focal point w_min=%.6f.\n', rep1.w_min);
    end
    fprintf('  %s\n\n', wmin_result.verdict);
end

fprintf('%s\n\n', bar);

end  

function print_instr_row(name, v0, v1, cfg)
    if ~any(strcmp(name, cfg.free))
        name = [name ' (FIXED)'];
    end
    fprintf('  %-28s  %14.6f  %14.6f\n', name, v0, v1);
end

function s = blocked_str(automation_share, eps_tol)
    if automation_share <= eps_tol
        s = 'Yes';
    else
        s = 'No';
    end
end
