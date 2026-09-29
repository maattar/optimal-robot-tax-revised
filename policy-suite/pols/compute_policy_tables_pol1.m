function [T1, T2, rep0, rep1] = compute_policy_tables_pol1(SPEC, x_opt, x_opt_label, x_bench_override)
% COMPUTE_POLICY_TABLES_POL1  
%   rep0, rep1  the underlying reports.
if nargin < 1 || isempty(SPEC), SPEC = 'baseline'; end
if nargin < 2, x_opt = []; end
if nargin < 3, x_opt_label = ''; end
if nargin < 4, x_bench_override = []; end
if ~isempty(x_opt) && isempty(x_opt_label)
    error(['compute_policy_tables_pol1: x_opt was given explicitly but x_opt_label was not']);
end

cfg = pol_config('pol1');
if ~cfg.built
    error('compute_policy_tables_pol1: pol1 is not marked built in pol_config.m.');
end
[cal_file, res_file, corner_file] = pol_filenames(cfg, SPEC);

% calib
if ~exist(cal_file, 'file')
    error('compute_policy_tables_pol1: %s not found.', cal_file);
end
S = load(cal_file);
zeta_star = S.zeta_star; fp = S.fp;
w_min_calib = zeta_star(5);

cp.alpha=zeta_star(1); cp.psi=zeta_star(2); cp.lam_bar=zeta_star(3); cp.eta_bar=zeta_star(4);
cp.delta=fp.delta; cp.theta=fp.theta; cp.nu=fp.nu; cp.sigma=fp.sigma;
cp.g_A=fp.g_A; cp.g_N=fp.g_N; cp.rho=fp.rho; cp.s_a1=fp.s_a1; cp.s_a2=fp.s_a2;
cp.s_n=1-fp.s_a1-fp.s_a2; cp.verbose=false;
cp.r     = fp.r;        % the capital tax never moves here, so no Euler step
cp.tau_K = fp.tau_K;    % fixed
cp.tau_L = fp.tau_L;    % fixed
cp.w_min = w_min_calib; % fixed





if isempty(x_opt)
    if exist(corner_file, 'file')
        Cres = load(corner_file);
        if ~isfield(Cres, 'corner')
            error(['compute_policy_tables_pol1: %s does not look like a ' ...
                   'locate_corner_pol1.m output.'], corner_file);
        end
        x_opt = Cres.corner;
        x_opt_label = sprintf('exact closed-form corner (from %s)', corner_file);
    elseif exist(res_file, 'file')
        R = load(res_file);
        x_opt = R.tau_R_star;
        x_opt_label = sprintf('SA optimum (from %s)', res_file);
    else
        error(['compute_policy_tables_pol1: neither %s nor %s found. Run ' ...
               'locate_corner_pol1(SPEC) (preferred -- exact) or ' ...
               'welfare_search_pol(''pol1'',SPEC) first, or pass x_opt ' ...
               'explicitly.'], corner_file, res_file);
    end
end

fprintf('\n===========================================================\n');
fprintf('  COMPUTE_POLICY_TABLES_POL1  SPEC=''%s''\n', SPEC);
fprintf('===========================================================\n');
fprintf('  Optimal-candidate point: tau_R=%.9f  (tau_K=%.6f, tau_L=%.6f, w_min=%.6f all FIXED)\n', ...
        x_opt, cp.tau_K, cp.tau_L, cp.w_min);
fprintf('  Label: %s\n\n', x_opt_label);

% The two reports: the status quo, at a zero robot tax with everything
% else at baseline, against the candidate
x_bench = 0.0;   
bench_label = 'calibrated status quo (tau_R=0)';
if ~isempty(x_bench_override)
    x_bench = x_bench_override;
    bench_label = 'x_bench_override';
end
rep0 = compute_policy_report_pol(x_bench, cfg, cp);
rep1 = compute_policy_report_pol(x_opt,   cfg, cp);
rep1.knife_edge_recovered = false;   


PLATEAU_EPS = 1e-6;
if ~rep1.converged
    rep1_nudged = compute_policy_report_pol(x_opt + PLATEAU_EPS, cfg, cp);
    if rep1_nudged.converged
        rep1_nudged.tau_R = x_opt;   % report the corner itself
        rep1_nudged.knife_edge_recovered = true;
        rep1 = rep1_nudged;
        fprintf('  NOTE: bare corner did not converge (regime=failed); recovered at corner+%.0e, same SWF (see rep1.knife_edge_recovered).\n', PLATEAU_EPS);
    end
end

if ~rep0.converged
    error(['compute_policy_tables_pol1: status-quo report (%s, tau_R=%.6f) ' ...
           'failed to converge (regime=%s). NOTE: the literal calibrated ' ...
           'benchmark (tau_R=0) sits in the AUTOMATION regime for this ' ...
           'calibration (w_min_calib > w_a1^{c,NA} at tau_R=0).'], bench_label, x_bench, rep0.regime);
end
if ~rep1.converged
    error(['compute_policy_tables_pol1: optimal-candidate report failed to ' ...
           'converge at tau_R=%.6f (regime=%s), even after retrying at ' ...
           'corner + %.0e. This is no longer the ordinary knife-edge case ' ...
           '(that is handled above) -- something else is wrong; inspect ' ...
           'solve_minwage_bgp/solve_joint_binding_bgp directly at this ' ...
           'tau_R before proceeding.'], x_opt, rep1.regime, PLATEAU_EPS);
end

% The two tables, through the shared printer every policy uses
[T1, T2, cev] = print_policy_tables_pol(rep0, rep1, cfg, cp, 'pol1', ...
    sprintf('candidate = %s', x_opt_label));

% Save
outfile = sprintf('policy_tables_result_%s_%s.mat', cfg.tag, SPEC);
save(outfile, 'SPEC', 'x_opt', 'x_opt_label', 'rep0', 'rep1', 'T1', 'T2', 'cev', 'cp');
fprintf('\nSaved %s\n\n', outfile);

end  
