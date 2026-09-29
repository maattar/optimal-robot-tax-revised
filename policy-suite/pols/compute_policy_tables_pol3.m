function [T1, T2, rep0, rep1] = compute_policy_tables_pol3(SPEC, x_opt, x_opt_label, x_bench_override)
% COMPUTE_POLICY_TABLES_POL3  The two comparison tables for Policy 3
%
%   [T1, T2, rep0, rep1] = compute_policy_tables_pol3(SPEC)
%   [T1, T2, rep0, rep1] = compute_policy_tables_pol3(SPEC, x_opt, label)
%
%   rep0, rep1  the underlying reports.


if nargin < 1 || isempty(SPEC), SPEC = 'baseline'; end
if nargin < 2, x_opt = []; end
if nargin < 3, x_opt_label = ''; end
if nargin < 4, x_bench_override = []; end
if ~isempty(x_opt) && isempty(x_opt_label)
    error(['compute_policy_tables_pol3: x_opt was given explicitly but x_opt_label was not']);
end

cfg = pol_config('pol3');
if ~cfg.built
    error('compute_policy_tables_pol3: pol3 is not marked built in pol_config.m.');
end
[cal_file, res_file, corner_file] = pol_filenames(cfg, SPEC);

% Load the calibration. Nothing here writes to it.
if ~exist(cal_file, 'file')
    error('compute_policy_tables_pol3: %s not found.', cal_file);
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
cp.tau_R = 0;           % pol3's defining feature (shared with pol2)
cp.w_min = w_min_calib; 






if isempty(x_opt)
    if exist(res_file, 'file')
        R = load(res_file);
        x_opt = [R.tau_K_star, R.tau_L_star];
        x_opt_label = sprintf('SA optimum (from %s)', res_file);
    elseif exist(corner_file, 'file')
        Cres = load(corner_file);
        if ~isfield(Cres, 'is_grid_quantized_only') || ~Cres.is_grid_quantized_only
            error(['compute_policy_tables_pol3: %s does not look like a locate_corner_pol3.m output.'], corner_file);
        end
        x_opt = [Cres.tau_K_star_grid, Cres.tau_L_star_grid];
        x_opt_label = sprintf(['GRID-QUANTIZED candidate (from %s) -NOT a verified optimum'], corner_file);
    else
        error(['compute_policy_tables_pol3: neither %s nor %s found. Run ' ...
               'welfare_search_pol(''pol3'',SPEC) or locate_corner_pol3(SPEC) ' ...
               'first, or pass x_opt explicitly.'], res_file, corner_file);
    end
end

fprintf('\n===========================================================\n');
fprintf('  COMPUTE_POLICY_TABLES_POL3  SPEC=''%s''\n', SPEC);
fprintf('===========================================================\n');
fprintf('  Optimal-candidate point: tau_K=%.6f  tau_L=%.6f  (tau_R=0, w_min=%.6f FIXED)\n', ...
        x_opt, cp.w_min);
fprintf('  Label: %s\n\n', x_opt_label);

% rep
x_bench = [fp.tau_K, fp.tau_L];
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
    x_opt_nudged = x_opt;
    x_opt_nudged(1) = x_opt_nudged(1) + PLATEAU_EPS;
    rep1_nudged = compute_policy_report_pol(x_opt_nudged, cfg, cp);
    if rep1_nudged.converged
        rep1_nudged.tau_K = x_opt(1);   % report the original candidate, not +eps
        rep1_nudged.knife_edge_recovered = true;
        rep1 = rep1_nudged;
        fprintf('      NOTE: bare candidate did not converge (regime=failed)', PLATEAU_EPS);
    end
end

if ~rep0.converged
    error(['compute_policy_tables_pol3: status-quo report (%s) failed to ' ...
           'converge at tau_K=%.6f tau_L=%.6f (regime=%s, w_min=%.6f FIXED).'], bench_label, x_bench, rep0.regime, cp.w_min);
end
if ~rep1.converged
    error(['compute_policy_tables_pol3: optimal-candidate report failed to converge at tau_K=%.6f tau_L=%.6f (regime=%s).'], ...
          x_opt, rep1.regime, PLATEAU_EPS);
end




[T1, T2, cev] = print_policy_tables_pol(rep0, rep1, cfg, cp, 'pol3', ...
    sprintf('candidate = %s', x_opt_label));

% Save
outfile = sprintf('policy_tables_result_%s_%s.mat', cfg.tag, SPEC);
save(outfile, 'SPEC', 'x_opt', 'x_opt_label', 'rep0', 'rep1', 'T1', 'T2', 'cev', 'cp');
fprintf('\nSaved %s\n\n', outfile);

end 
