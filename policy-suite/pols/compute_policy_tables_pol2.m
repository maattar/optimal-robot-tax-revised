function [T1, T2, rep0, rep1] = compute_policy_tables_pol2(SPEC, x_opt, x_opt_label, x_bench_override)
% COMPUTE_POLICY_TABLES_POL2  The two comparison tables for pol2
%
%   [T1, T2, rep0, rep1] = compute_policy_tables_pol2(SPEC)
%   [T1, T2, rep0, rep1] = compute_policy_tables_pol2(SPEC, x_opt, label)
%
%   rep0, rep1  the underlying reports.
if nargin < 1 || isempty(SPEC), SPEC = 'baseline'; end
if nargin < 2, x_opt = []; end
if nargin < 3, x_opt_label = ''; end
if nargin < 4, x_bench_override = []; end
if ~isempty(x_opt) && isempty(x_opt_label)
    error(['compute_policy_tables_pol2: x_opt was given explicitly but ' ...
           'x_opt_label was not.']);
end

cfg = pol_config('pol2');
if ~cfg.built
    error('compute_policy_tables_pol2: pol2 is not marked built in pol_config.m.');
end
[cal_file, res_file, corner_file] = pol_filenames(cfg, SPEC);

% Load the calinration
if ~exist(cal_file, 'file')
    error('compute_policy_tables_pol2: %s not found.', cal_file);
end
S = load(cal_file);
zeta_star = S.zeta_star; fp = S.fp;
w_min_calib = zeta_star(5);

cp.alpha=zeta_star(1); cp.psi=zeta_star(2); cp.lam_bar=zeta_star(3); cp.eta_bar=zeta_star(4);
cp.delta=fp.delta; cp.theta=fp.theta; cp.nu=fp.nu; cp.sigma=fp.sigma;
cp.g_A=fp.g_A; cp.g_N=fp.g_N; cp.rho=fp.rho; cp.s_a1=fp.s_a1; cp.s_a2=fp.s_a2;
cp.s_n=1-fp.s_a1-fp.s_a2; cp.verbose=false;
[Phi,~] = phi_from_calibration(fp,false);
cp.Phi = Phi;
cp.tau_R = 0;           
cp.tau_L = fp.tau_L;    
cp.w_min = w_min_calib; 




if isempty(x_opt)
    if exist(res_file, 'file')
        R = load(res_file);
        x_opt = R.tau_K_star;
        x_opt_label = sprintf('SA optimum (from %s)', res_file);
    elseif exist(corner_file, 'file')
        Cres = load(corner_file);
        if ~isfield(Cres, 'is_grid_quantized_only') || ~Cres.is_grid_quantized_only
            error(['compute_policy_tables_pol2: %s does not look like a ' ...
                   'locate_corner_pol2.m output.'], corner_file);
        end
        x_opt = Cres.tau_K_star_grid;
        x_opt_label = sprintf(['GRID-QUANTIZED candidate (from %s) -- ' ...
                                'NOT a verified optimum, see that file''s header'], corner_file);
    else
        error(['compute_policy_tables_pol2: neither %s nor %s found. Run ' ...
               'welfare_search_pol(''pol2'',SPEC) or locate_corner_pol2(SPEC) ' ...
               'first, or pass x_opt explicitly.'], res_file, corner_file);
    end
end

fprintf('\n===========================================================\n');
fprintf('  COMPUTE_POLICY_TABLES_POL2  SPEC=''%s''\n', SPEC);
fprintf('===========================================================\n');
fprintf('  Optimal-candidate point: tau_K=%.6f  (tau_R=0, tau_L=%.6f, w_min=%.6f FIXED)\n', ...
        x_opt, cp.tau_L, cp.w_min);
fprintf('  Label: %s\n\n', x_opt_label);

% The two reports: the status quo, at the calibrated capital tax with no
% robot tax and the other two fixed, against the candidate
x_bench = fp.tau_K;
bench_label = 'calibrated status quo';
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
        rep1_nudged.tau_K = x_opt;   % report the original candidate, not +eps
        rep1_nudged.knife_edge_recovered = true;
        rep1 = rep1_nudged;
        fprintf('  NOTE: bare candidate did not converge (regime=failed); recovered at tau_K+%.0e, same SWF (see rep1.knife_edge_recovered).\n', PLATEAU_EPS);
    end
end

if ~rep0.converged
    error(['compute_policy_tables_pol2: status-quo report (%s, tau_K=%.6f) ' ...
           'failed to converge (regime=%s). This calibration''s status quo ' ...
           'falls in the automation regime at tau_R=0'], bench_label, x_bench, rep0.regime);
end
if ~rep1.converged
    error(['compute_policy_tables_pol2: optimal-candidate report failed to ' ...
           'converge at tau_K=%.6f (regime=%s).'], x_opt, rep1.regime, PLATEAU_EPS);
end

% The two tables
[T1, T2, cev] = print_policy_tables_pol(rep0, rep1, cfg, cp, 'pol2', ...
    sprintf('candidate = %s', x_opt_label));

% Save
outfile = sprintf('policy_tables_result_%s_%s.mat', cfg.tag, SPEC);
save(outfile, 'SPEC', 'x_opt', 'x_opt_label', 'rep0', 'rep1', 'T1', 'T2', 'cev', 'cp');
fprintf('\nSaved %s\n\n', outfile);

end 