% CALIBRATE_PARALLEL  Outer calibration loop.
%
%   Script; needs the Global Optimization and Parallel Computing toolboxes.
%   Phase 1 evaluates N_SCAN random parameter vectors in parallel, Phase 2
%   runs simulated annealing from the best N_REFINE of them, one run per
%   worker.
%
%   IMPORTANT: The algorithm does not find *the* exact minimum in each run
%              because refinement is stochastic. See Online Appendix D for
%              details.
%
%   Free parameters   zeta = [alpha, psi, lam_bar, eta_bar, w_min]
%
%   Targeted moment            Source        
%   -----------------------------------------
%   (K+R)/Y                    PWT 10.01     
%   Labour income share        PWT 10.01     
%   Robot investment share     ACES 2022     
%   P50/P10 wage ratio         BLS Table 5   
%   Aggregate unemployment     CBO      
%
%   The targets themselves are set in M_data below, along with the
%   validation moments that are computed and reported but not matched.
%   Five moments for five parameters: the system is exactly identified.
%
%   w_min is read as an effective economy-wide wage floor rather than the
%   statutory minimum. The OECD min-to-mean wage ratio is untargeted, and 
%   the model is expected to come out above the OECD figure for the 
%   statutory minimum. Rather than impose that as a hard constraint, the 
%   objective carries a soft penalty for falling below it.
%
%   The benchmark requires the floor to bind for type-1 workers, and
%   solve_bgp enforces that for every draw: the competitive BGP must solve
%   and confirm its regime, w_min must strictly exceed the competitive
%   type-1 wage, and the min-wage BGP must converge. Any failure returns 
%   bgp_ok = false and the draw is penalised in calibration_obj.
%
%   Each SPEC writes calibration_result_<SPEC>.mat, so runs never collide.

clearvars

% =========================================================================
% SPECIFICATION
% =========================================================================
%   The robustness alternatives differ from the baseline only in the
%   fixed parameters set in the switch below.  Targets, bounds, algorithm
%   settings and the objective are identical across all four. 
% 
%   Run the script once per SPEC.
%
%     'baseline'      benchmark
%     'auto33'        larger automatable population, s_a1 + s_a2 = 1/3,
%                     with only s_a2 moving
%     'equal_shares'  same automatable mass as the benchmark, split evenly
%                     between the two types
%     'acemoglu2020'  income tax rates from Acemoglu, Manera and Restrepo
%                     (2020)
% =========================================================================
SPEC = 'baseline';

VALID_SPECS = {'baseline', 'auto33', 'equal_shares', 'acemoglu2020'};
if ~ismember(SPEC, VALID_SPECS)
    error('calibrate_parallel: unknown SPEC ''%s''. Valid options: %s.', ...
          SPEC, strjoin(VALID_SPECS, ', '));
end

% =========================================================================
% FIXED PARAMETERS
% =========================================================================
fp.r      = 0.07;
fp.delta  = 0.04317;
fp.theta  = 1 - (1/0.704);
fp.tau_K  = 0.077180756;
fp.tau_L  = 0.187149534;
fp.nu     = 0.46;
fp.g_A    = 0.01141;
fp.g_N    = 0.005859;
fp.rho    = 0.030;

fp.s_a1   = 0.10;
fp.s_a2   = 0.18;

switch SPEC
    case 'baseline'
        % no overrides

    case 'auto33'
        fp.s_a1 = 0.10;          % unchanged
        fp.s_a2 = 1/3 - 0.10;    % s_a1 + s_a2 = 1/3 exactly

    case 'equal_shares'
        fp.s_a1 = 0.14;
        fp.s_a2 = 0.14;          % same total as the benchmark

    case 'acemoglu2020'
        fp.tau_K = 0.100;
        fp.tau_L = 0.255;
end

fp.s_n = 1 - fp.s_a1 - fp.s_a2;

% Each spec has its own expected s_n.
switch SPEC
    case {'baseline', 'acemoglu2020', 'equal_shares'}
        s_n_expected = 0.72;
    case 'auto33'
        s_n_expected = 1 - 1/3;
end
if abs(fp.s_n - s_n_expected) > 1e-9
    error(['calibrate: fp.s_n = %.6f does not match the expected value ' ...
           '%.6f for SPEC=''%s''; check fp.s_a1/fp.s_a2 overrides.'], ...
          fp.s_n, s_n_expected, SPEC);
end

fprintf('\n>>> calibrate_parallel running under SPEC = ''%s'' <<<\n', SPEC);
fprintf('    fp.nu=%.4f  fp.tau_K=%.6f  fp.tau_L=%.6f  fp.s_a1=%.4f  fp.s_a2=%.4f\n\n', ...
        fp.nu, fp.tau_K, fp.tau_L, fp.s_a1, fp.s_a2);

% Set true to identify eta_bar off the bottom-quintile income share
% instead of the P50/P10 wage ratio.
USE_BQ_FOR_ETA = false;

% sigma follows from the Euler equation on the BGP.
fp.sigma = (fp.r * (1 - fp.tau_K) - fp.g_N - fp.rho) / fp.g_A;
fprintf('Derived sigma = %.4f\n', fp.sigma);
if fp.sigma <= 0
    error('calibrate: sigma <= 0; check r, tau_K, g_N, rho, g_A.');
end

% =========================================================================
% DATA TARGETS
% =========================================================================

M_data.K2Y            = 3.498;   % (K+R)/Y, PWT 10.01                      [targeted]
M_data.labor_share    = 0.594;   % pre-tax labour income share, PWT 10.01  [targeted]
M_data.IRshare        = 0.011;   % robot investment share, ACES 2022       [targeted]
M_data.p50_p10        = 2.09;    % BLS Usual Weekly Earnings, Table 5      [targeted by default]
M_data.bq_share       = 0.075;   % lowest quintile, after T&T, CBO 58353 Exh. 20
                                 %                                         [targeted if USE_BQ_FOR_ETA]
M_data.u_a1_nairu     = 0.044;   % CBO Pub. 55551                          [targeted]

% --- Validation only ---
M_data.min_mean_ratio = 0.256772;   % OECD MIN2AVE
M_data.gini_before    = 0.517;      % CBO Pub. 58353; their "before T&T"
                                    % includes social insurance rather than
                                    % being pure market income
M_data.gini_after     = 0.4323;     % CBO Pub. 58353
M_data.bq_share_before    = 0.0378;   % lowest quintile, before T&T
M_data.top20_share_before = 0.5470;   % highest quintile, before T&T
M_data.top20_share_after  = 0.4834;   % highest quintile, after T&T

% =========================================================================
% BOUNDS AND ALGORITHM SETTINGS
% =========================================================================
%   zeta = [alpha, psi, lam_bar, eta_bar, w_min]
lb = [0.01, 0.01, 0.01, 1.01, 0.01];
ub = [0.50, 0.25, 20.0, 20.00, 2.00];

N_SCAN   = 100000;
N_REFINE =    200;

sa_opts = optimoptions('simulannealbnd', ...
    'Display',               'off',           ...
    'MaxIterations',          5000,            ...
    'MaxFunctionEvaluations', 20000,           ...
    'InitialTemperature',     0.1,             ...
    'ReannealInterval',       300,             ...
    'TemperatureFcn',        'temperatureexp', ...
    'AnnealingFcn',          'annealingfast',  ...
    'ObjectiveLimit',         1e-8);

pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local');
end
fprintf('Parallel pool: %d workers.\n\n', pool.NumWorkers);

% =========================================================================
% PHASE 1 --- COARSE SCAN
% =========================================================================
fprintf('=== Phase 1: Coarse scan (%d evaluations, parallel) ===\n', N_SCAN);

rng(42, 'twister');
zeta_draws = zeros(N_SCAN, 5);
for i = 1:N_SCAN
    zeta_draws(i,:) = draw_feasible(lb, ub);
end

Q_scan = inf(N_SCAN, 1);
parfor i = 1:N_SCAN
    Q_scan(i) = calibration_obj(zeta_draws(i,:), fp, M_data, USE_BQ_FOR_ETA);  %#ok<PFBNS>
end

n_fail = sum(~isfinite(Q_scan));
[Q_sorted, sort_idx] = sort(Q_scan);

fprintf('Phase 1 complete.  Failures: %d / %d\n', n_fail, N_SCAN);
fprintf('Best Q = %.6e  at  alpha=%.4f  psi=%.4f  lam_bar=%.4f  eta_bar=%.4f  w_min=%.4f\n\n', ...
        Q_sorted(1), zeta_draws(sort_idx(1),:));

% =========================================================================
% PHASE 2 --- SA REFINEMENT
% =========================================================================
fprintf('=== Phase 2: SA refinement (%d runs, parallel) ===\n', N_REFINE);

n_candidates = min(N_REFINE, sum(isfinite(Q_sorted)));
z0_all = zeta_draws(sort_idx(1:n_candidates), :);
Q0_all = Q_scan(sort_idx(1:n_candidates));

Q_ref    = inf(n_candidates, 1);
zeta_ref = zeros(n_candidates, 5);

parfor j = 1:n_candidates
    
    z0_j = z0_all(j,:);   %#ok<PFBNS>
    Q0_j = Q0_all(j);     %#ok<PFBNS>
    try
        [z_opt, Q_opt] = simulannealbnd( ...
            @(z) calibration_obj(z, fp, M_data, USE_BQ_FOR_ETA), ...
            z0_j, lb, ub, sa_opts);               %#ok<PFBNS>
        Q_ref(j)      = Q_opt;
        zeta_ref(j,:) = z_opt;
    catch
        % Keep the Phase-1 value rather than losing the start point.
        Q_ref(j)      = Q0_j;
        zeta_ref(j,:) = z0_j;
    end
end

fprintf('Phase 2 results:\n');
for j = 1:n_candidates
    fprintf('  Run %2d: Q = %.6e  (alpha=%.4f  psi=%.4f  eta_bar=%.4f  w_min=%.4f)\n', ...
            j, Q_ref(j), zeta_ref(j,1), zeta_ref(j,2), zeta_ref(j,4), zeta_ref(j,5));
end

% =========================================================================
% BEST RUN
% =========================================================================
[Q_star, best_j] = min(Q_ref);
zeta_star        = zeta_ref(best_j,:);

fprintf('\n=== Calibration complete ===\n');
fprintf('Global minimum: Q* = %.6e\n', Q_star);
fprintf('zeta*:  alpha=%.6f  psi=%.6f  lam_bar=%.6f  eta_bar=%.6f  w_min=%.6f\n', ...
        zeta_star);
fprintf('Fixed:  s_a1=%.4f  s_a2=%.4f  s_n=%.4f  (not searched)\n\n', ...
        fp.s_a1, fp.s_a2, fp.s_n);

% Re-solve at the optimum to recover the full equilibrium.
params_star = build_params(zeta_star, fp);
[sol_comp_star, sol_mw_star, bgp_ok_star] = solve_bgp(params_star);

if bgp_ok_star
    sc_star   = compute_scalars(params_star);
    moms_star = compute_moments(sol_comp_star, sol_mw_star, sc_star, params_star);
    display_results(zeta_star, fp, moms_star, M_data, USE_BQ_FOR_ETA, SPEC);
else
    warning('calibrate: BGP did not converge on final re-evaluation at zeta*.');
end

outfile = sprintf('calibration_result_%s.mat', SPEC);
save(outfile, ...
     'zeta_star', 'Q_star', 'fp', 'M_data', 'USE_BQ_FOR_ETA', ...
     'sol_comp_star', 'sol_mw_star', 'moms_star', ...
     'Q_scan', 'zeta_draws', 'Q_ref', 'zeta_ref', 'SPEC');
fprintf('Results saved to %s\n', outfile);


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function Q = calibration_obj(zeta, fp, M_data, use_bq_for_eta)
% Sum of squared relative deviations across the five targeted moments,
% plus the soft floor penalty.

    PENALTY = 1e10;

    % alpha + psi < 1 is needed for the exponent on s_n in A.
    if (zeta(1) + zeta(2)) >= 0.99
        Q = PENALTY;
        return
    end

    params = build_params(zeta, fp);

    [sol_comp, sol_mw, bgp_ok] = solve_bgp(params);
    if ~bgp_ok
        Q = PENALTY;
        return
    end

    sc   = compute_scalars(params);
    moms = compute_moments(sol_comp, sol_mw, sc, params);

    if ~isfinite(moms.min_mean_ratio)
        Q = PENALTY;
        return
    end

    if use_bq_for_eta
        eta_model = moms.bq_share_after;
        eta_data  = M_data.bq_share;
    else
        eta_model = moms.p50_p10;
        eta_data  = M_data.p50_p10;
    end

    m_vec = [moms.K2Y;   moms.labor_share;   moms.IRshare; ...
             eta_model;  moms.u_a1];
    d_vec = [M_data.K2Y; M_data.labor_share; M_data.IRshare; ...
             eta_data;   M_data.u_a1_nairu];

    if any(~isfinite(m_vec))
        Q = PENALTY;
        return
    end

    Q_fit = sum( ((d_vec - m_vec) ./ d_vec).^2 );

    PEN_WEIGHT = 50;
    violation  = max(0, (M_data.min_mean_ratio - moms.min_mean_ratio) / M_data.min_mean_ratio);
    Q_penalty  = PEN_WEIGHT * violation^2;

    Q = Q_fit + Q_penalty;
end


function params = build_params(zeta, fp)
    params.alpha   = zeta(1);
    params.psi     = zeta(2);
    params.lam_bar = zeta(3);
    params.eta_bar = zeta(4);
    params.w_min   = zeta(5);

    params.s_a1    = fp.s_a1;
    params.s_a2    = fp.s_a2;
    params.s_n     = fp.s_n;

    params.r       = fp.r;
    params.delta   = fp.delta;
    params.theta   = fp.theta;
    params.tau_K   = fp.tau_K;
    params.tau_L   = fp.tau_L;
    params.nu      = fp.nu;
    params.sigma   = fp.sigma;
    params.g_A     = fp.g_A;
    params.g_N     = fp.g_N;
    params.rho     = fp.rho;

    params.tau_R   = 0;
    params.verbose = false;
end


function zeta = draw_feasible(lb, ub)
% Uniform draw on [lb, ub], rejecting anything with alpha + psi >= 0.99.
    while true
        zeta = lb + (ub - lb) .* rand(1, 5);
        if (zeta(1) + zeta(2)) < 0.99
            return
        end
    end
end


function display_results(zeta_star, fp, moms, M_data, use_bq_for_eta, spec)

    fprintf('\n====================================================\n');
    fprintf('  CALIBRATION RESULTS \n');
    fprintf('  SPEC = ''%s''\n', spec);
    fprintf('  w_min = effective wage floor\n');
    fprintf('====================================================\n');

    fprintf('\nFixed parameters:\n');
    fprintf('  r      = %.4f   delta  = %.4f   theta  = %.4f\n', fp.r, fp.delta, fp.theta);
    fprintf('  tau_K  = %.4f   tau_L  = %.4f   nu     = %.4f\n', fp.tau_K, fp.tau_L, fp.nu);
    fprintf('  g_A    = %.4f   g_N    = %.4f   rho    = %.4f\n', fp.g_A, fp.g_N, fp.rho);
    fprintf('  sigma  = %.4f\n', fp.sigma);
    fprintf('  s_a1   = %.4f   s_a2   = %.4f   s_n    = %.4f   (fixed)\n', ...
            fp.s_a1, fp.s_a2, fp.s_n);

    fprintf('\nCalibrated parameters (zeta*):\n');
    fprintf('  alpha   = %.6f\n', zeta_star(1));
    fprintf('  psi     = %.6f\n', zeta_star(2));
    fprintf('  lam_bar = %.6f\n', zeta_star(3));
    fprintf('  eta_bar = %.6f\n', zeta_star(4));
    fprintf('  w_min   = %.6f\n', zeta_star(5));

    fprintf('\nTargeted moments (5, exactly identified):\n');
    fprintf('  %-35s  %8s  %8s  %8s\n', 'Moment', 'Data', 'Model', '% Error');
    fprintf('  %s\n', repmat('-', 1, 65));

    if use_bq_for_eta
        eta_name  = 'Bottom-quintile share (post-T&T)';
        eta_model = moms.bq_share_after;
        eta_data  = M_data.bq_share;
    else
        eta_name  = 'P50/P10 wage ratio';
        eta_model = moms.p50_p10;
        eta_data  = M_data.p50_p10;
    end

    tgt_names = {'(K+R)/Y',                   ...
                 'Labour income share',        ...
                 'Robot inv. share',           ...
                 eta_name,                     ...
                 'Aggregate unemployment'};
    tgt_model = [moms.K2Y;          moms.labor_share;    moms.IRshare; ...
                 eta_model;         moms.u_a1];
    tgt_data  = [M_data.K2Y;        M_data.labor_share;  M_data.IRshare; ...
                 eta_data;          M_data.u_a1_nairu];

    for i = 1:5
        err = 100 * (tgt_model(i) - tgt_data(i)) / tgt_data(i);
        fprintf('  %-35s  %8.4f  %8.4f  %7.2f%%\n', ...
                tgt_names{i}, tgt_data(i), tgt_model(i), err);
    end

    fprintf('\nValidation moments (untargeted):\n');

    fprintf('\n  Distribution -- quintile income shares and Gini:\n');
    fprintf('  %-36s  %10s  %10s\n', 'Moment', 'Data', 'Model');
    fprintf('  %s\n', repmat('-', 1, 61));

    row('Bottom quintile share, before T&T', M_data.bq_share_before,    moms.bq_share_before);
    if use_bq_for_eta
        fprintf('  %-36s  %10s  %10s\n', 'Bottom quintile share, after T&T', '[targeted]', '[see above]');
    else
        row('Bottom quintile share, after T&T',  M_data.bq_share,       moms.bq_share_after);
    end
    row('Top quintile share, before T&T',    M_data.top20_share_before, moms.top20_share_before);
    row('Top quintile share, after T&T',     M_data.top20_share_after,  moms.top20_share_after);
    fprintf('  %s\n', repmat('-', 1, 61));
    row('Gini before T&T',                   M_data.gini_before,        moms.gini_before);
    row('Gini after T&T',                    M_data.gini_after,         moms.gini_after);
    row('Gini change (after - before)',      M_data.gini_after - M_data.gini_before, ...
                                             moms.gini_after - moms.gini_before);
    fprintf('      (negative = redistribution compresses the distribution)\n');

    fprintf('\n  Wages and employment:\n');
    fprintf('  %-36s  %10s  %10s\n', 'Moment', 'Data', 'Model');
    fprintf('  %s\n', repmat('-', 1, 61));
    row('Min-to-mean wage ratio',            M_data.min_mean_ratio,     moms.min_mean_ratio);
    fprintf('      (expect model ABOVE %.3f -- effective floor > federal min.)\n', M_data.min_mean_ratio);
    fprintf('  %-36s  %10s  %10.4f\n', 'Wage ratio w_a2/w_a1',          '>1 (no pt.)', moms.wage_ratio);
    fprintf('  %-36s  %10s  %10.4f\n', 'Employed type-1 share (=L_a1)', 'no data cp.', moms.emp_minwage_sh);
    if use_bq_for_eta
        row('P50/P10 wage ratio',            M_data.p50_p10,            moms.p50_p10);
    else
        fprintf('  %-36s  %10s  %10s\n', 'P50/P10 wage ratio', '[targeted]', '[see above]');
    end

    fprintf('\nSocial welfare: SWF = %.6f\n', moms.SWF);
    fprintf('====================================================\n\n');
end


function row(name, data_val, model_val)

    if isnan(data_val)
        fprintf('  %-36s  %10s  %10.4f\n', name, 'pending', model_val);
    else
        fprintf('  %-36s  %10.4f  %10.4f\n', name, data_val, model_val);
    end
end
