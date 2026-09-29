function [T1, T2, rep0, rep1] = compute_policy_tables_pol5(SPEC, x_opt, x_opt_label, wmin_result)
% COMPUTE_POLICY_TABLES_POL5  The two comparison tables for pol5.

%   [T1, T2, rep0, rep1] = compute_policy_tables_pol5(SPEC)
%   [T1, T2, rep0, rep1] = compute_policy_tables_pol5(SPEC, x_opt, label)
%


if nargin < 1 || isempty(SPEC), SPEC = 'baseline'; end
if nargin < 2, x_opt = []; end
if nargin < 3, x_opt_label = ''; end
if nargin < 4, wmin_result = []; end
if ~isempty(x_opt) && isempty(x_opt_label)
    error(['compute_policy_tables_pol5: x_opt was given explicitly but ' ...
           'x_opt_label was not -- pass a short description of what this ' ...
           'point is.']);
end

cfg = pol_config('pol5');
if ~cfg.built
    error('compute_policy_tables_pol5: pol5 is not marked built in pol_config.m.');
end
[cal_file, ~, corner_file] = pol_filenames(cfg, SPEC);

% Load the calibration. Nothing here writes to it.
if ~exist(cal_file, 'file')
    error('compute_policy_tables_pol5: %s not found.', cal_file);
end
S = load(cal_file);
zeta_star = S.zeta_star; fp = S.fp;
w_min_calib = zeta_star(5);

cp.alpha=zeta_star(1); cp.psi=zeta_star(2); cp.lam_bar=zeta_star(3); cp.eta_bar=zeta_star(4);
cp.delta=fp.delta; cp.theta=fp.theta; cp.nu=fp.nu; cp.sigma=fp.sigma;
cp.g_A=fp.g_A; cp.g_N=fp.g_N; cp.rho=fp.rho; cp.s_a1=fp.s_a1; cp.s_a2=fp.s_a2;
cp.s_n=1-fp.s_a1-fp.s_a2; cp.verbose=false;
cp.r     = fp.r;        
cp.tau_K = fp.tau_K;    
cp.tau_L = fp.tau_L;    
cp.tau_R = 0;           

if isempty(x_opt)
    if ~exist(corner_file, 'file')
        error(['compute_policy_tables_pol5: %s not found. Run ' ...
               'locate_corner_pol5(SPEC) first, or pass x_opt ' ...
               'explicitly.'], corner_file);
    end
    Cres = load(corner_file);
    if ~isfield(Cres, 'w_min_star')
        error(['compute_policy_tables_pol5: %s does not look like a ' ...
               'locate_corner_pol5.m output.'], corner_file);
    end
    x_opt = Cres.w_min_star;
    x_opt_label = sprintf('focal point of the flat interval (from %s)', corner_file);
end


if ~isempty(wmin_result) && isfield(wmin_result, 'is_indeterminate') ...
        && wmin_result.is_indeterminate && x_opt ~= 0
    fprintf(['  NOTE: wmin_identification returned INDETERMINATE; ' ...
             'reporting the floor at its focal value 0 rather than ' ...
             '%.9f (welfare is identical).\n'], x_opt);
    x_opt = 0;
    x_opt_label = [x_opt_label ' -> forced to focal 0 (undetermined)'];
end

fprintf('\n===========================================================\n');
fprintf('  COMPUTE_POLICY_TABLES_POL5  SPEC=''%s''\n', SPEC);
fprintf('===========================================================\n');
fprintf('  Optimal-candidate point: w_min=%.9f  (tau_R=%.6f, tau_K=%.6f, tau_L=%.6f all FIXED)\n', ...
        x_opt, cp.tau_R, cp.tau_K, cp.tau_L);
fprintf('  Label: %s\n\n', x_opt_label);


x_bench = w_min_calib;
bench_label = 'calibrated status quo (w_min = w_min_calib, tau_R = 0)';

rep0 = compute_policy_report_pol(x_bench, cfg, cp);
rep1 = compute_policy_report_pol(x_opt,   cfg, cp);

if ~rep0.converged
    error(['compute_policy_tables_pol5: status-quo report (%s, ' ...
           'w_min=%.6f) failed to converge (regime=%s). The calibrated ' ...
           'benchmark sits in the AUTOMATION regime with the floor ' ...
           'binding, which needs fsolve+optimoptions (MATLAB-only).'], ...
          bench_label, x_bench, rep0.regime);
end
if ~rep1.converged
    error(['compute_policy_tables_pol5: optimal-candidate report failed ' ...
           'to converge at w_min=%.6f (regime=%s). A floor at or below ' ...
           'the threshold cannot bind, so the competitive solver should ' ...
           'settle it; inspect solve_competitive_bgp directly at this ' ...
           'floor before proceeding.'], x_opt, rep1.regime);
end


[T1, T2, cev] = print_policy_tables_pol(rep0, rep1, cfg, cp, 'pol5', ...
    sprintf('candidate = %s', x_opt_label));

% Save
outfile = sprintf('policy_tables_result_%s_%s.mat', cfg.tag, SPEC);
save(outfile, 'SPEC', 'x_opt', 'x_opt_label', 'rep0', 'rep1', 'T1', 'T2', 'cev', 'cp');
fprintf('\nSaved %s\n\n', outfile);

end  
