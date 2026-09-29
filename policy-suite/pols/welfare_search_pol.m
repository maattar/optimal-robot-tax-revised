function welfare_search_pol(pol_id, SPEC, bounds, N_GRID, N_REFINE)
% WELFARE_SEARCH_POL  Searches for the optimal policy
%
%   welfare_search_pol('pol6')
%   welfare_search_pol('pol6', 'baseline', bounds)
%   welfare_search_pol('pol4', 'baseline', bounds, [50 50 50], 100)

if nargin < 2 || isempty(SPEC),   SPEC   = 'baseline'; end
if nargin < 3,                    bounds = struct();   end
if nargin < 4,                    N_GRID = [];          end
if nargin < 5,                    N_REFINE = [];        end

cfg = pol_config(pol_id);
if ~cfg.built
    error(['welfare_search_pol: %s is not yet built   '], pol_id);
end
[cal_file, res_file] = pol_filenames(cfg, SPEC);

% CALIBRATION 
if ~exist(cal_file, 'file')
    error(['welfare_search_pol: %s not found. Run calibrate_parallel.m ' ...
           'first, or copy the .mat file into this folder.'], cal_file);
end
S = load(cal_file);
if ~isfield(S, 'zeta_star') || ~isfield(S, 'fp')
    error('welfare_search_pol: %s must contain both ''zeta_star'' and ''fp''.', cal_file);
end
if isfield(S, 'SPEC') && ~isempty(SPEC) && ~strcmp(S.SPEC, SPEC)
    error(['welfare_search_pol: SPEC MISMATCH. Requested ''%s'' but %s ' ...
           'identifies itself as ''%s''.'], SPEC, cal_file, S.SPEC);
end
zeta_star = S.zeta_star;  fp = S.fp;
w_min_calib = zeta_star(5);

fprintf('\n===========================================================\n');
fprintf('  WELFARE_SEARCH_POL(''%s'')  SPEC=''%s''  free=%s\n', ...
        pol_id, SPEC, strjoin(cfg.free, ', '));
fprintf('===========================================================\n');
fprintf('  alpha=%.6f  psi=%.6f  lam_bar=%.6f  eta_bar=%.6f  w_min(calib)=%.6f\n', ...
        zeta_star);

% =========================================================================
% The parameters
% =========================================================================
cp.alpha   = zeta_star(1);
cp.psi     = zeta_star(2);
cp.lam_bar = zeta_star(3);
cp.eta_bar = zeta_star(4);
cp.delta = fp.delta;  cp.theta = fp.theta;  cp.nu   = fp.nu;
cp.sigma = fp.sigma;  cp.g_A   = fp.g_A;    cp.g_N  = fp.g_N;
cp.rho   = fp.rho;    cp.s_a1  = fp.s_a1;   cp.s_a2 = fp.s_a2;
cp.s_n   = 1 - fp.s_a1 - fp.s_a2;
cp.verbose = false;

if cfg.needs_euler
    [Phi, ~] = phi_from_calibration(fp, true);
    cp.Phi = Phi;
    r_base_recon = euler_closure(Phi, fp.tau_K);
    nest_gap = abs(r_base_recon - fp.r);
    % Phi is formed as fp.r*(1-fp.tau_K) and the rate recovered by dividing
    % the same factor back out, so the round trip is exact only where that
    % multiply and divide happen to cancel in double precision. They cancel
    % at the calibrated tau_K and not at the round 0.1 the acemoglu2020
    % spec sets, where the two land one ulp apart. Demanding bit-equality
    % here therefore tests the arithmetic rather than the calibration. The
    % tolerance is locate_corner_pol's, which checks the same quantity.
    NEST_TOL = 1e-9;
    if nest_gap > NEST_TOL
        error(['welfare_search_pol: NESTING GUARD FAILED. r(tau_K_base) ' ...
               '= %.17g does not reproduce fp.r = %.17g (gap %.3e, ' ...
               'tolerance %.1e). Refusing to run a search whose baseline ' ...
               'does not nest the completed results.'], ...
              r_base_recon, fp.r, nest_gap, NEST_TOL);
    end
    if nest_gap == 0
        fprintf('  Euler round-trip: r(tau_K_base) == fp.r bit-exact.  \n');
    else
        fprintf(['  Euler round-trip: r(tau_K_base) reproduces fp.r to ' ...
                 '%.3e (%.1f ulp).  []\n'], ...
                nest_gap, nest_gap / eps(fp.r));
    end
else
    cp.r = fp.r;
end

% Fixed (non-free) policy instruments, at their baseline/calibrated values.
if ~any(strcmp('tau_K', cfg.free)), cp.tau_K = fp.tau_K; end
if ~any(strcmp('tau_L', cfg.free)), cp.tau_L = fp.tau_L; end
if ~any(strcmp('w_min', cfg.free)), cp.w_min = w_min_calib; end
if ~any(strcmp('tau_R', cfg.free)), cp.tau_R = 0; end
% Fixed at tau_R=0 (not the calibrated fp.tau_R, which does not exist --
% tau_R is a POLICY instrument, never part of the structural calibration).
% This matters for any pol where tau_R is fixed rather than free (pol2,
% pol3): without this line, eval_policy_pol would receive a cp with no
% .tau_R field at all.

% =========================================================================
% BOUNDS  (defaults reproduce each original file's documented box)
% =========================================================================
DEFAULT_BOUNDS.tau_R = [-0.10, 3.0];
DEFAULT_BOUNDS.tau_K = [-1.00, 1.0];
DEFAULT_BOUNDS.tau_L = [-0.50, 1.0];
DEFAULT_BOUNDS.w_min = [0.0, 3.0 * w_min_calib];

nfree = numel(cfg.free);
lb = zeros(1, nfree); ub = zeros(1, nfree);
for i = 1:nfree
    fn = cfg.free{i};
    if isfield(bounds, fn)
        b = bounds.(fn);
    else
        b = DEFAULT_BOUNDS.(fn);
    end
    lb(i) = b(1); ub(i) = b(2);
    fprintf('  bound %-6s : [%.4f, %.4f]%s\n', fn, b(1), b(2), ...
            tern(isfield(bounds, fn), '', '  (default)'));
end

if isempty(N_GRID)
    % Default grid budgets, tuned per pol so that Phase 1 stays within a
    % reasonable wall-clock target while giving each pol's search space
    % enough resolution to reliably bracket the optimum before Phase 2
    % (simulated annealing) refines it. The multi-D pols share a common
    % TARGET_TOTAL_GRID (total grid points, split evenly per dimension);
    % pol1/pol2 are handled separately since both are 1-D.
    if strcmp(pol_id, 'pol2')
        % pol2 is 1-D (tau_K only) like pol1, but unlike pol1 there is no
        % closed-form optimum and no flat plateau -- SWF genuinely varies
        % with tau_K on both sides of tau_K_flip (see locate_corner_pol2.m),
        % so this grid needs finer resolution than pol1's plateau
        % cross-check to usefully bracket the true optimum before Phase 2
        % SA refines it.
        N_GRID = 5000;
    elseif strcmp(pol_id, 'pol1')
        % pol1 is 1-D (tau_R only, everything else fixed including w_min
        % -- unlike pol6's 2-D case). The exact optimum is already known
        % in closed form (locate_corner_pol1.m); this grid is purely a
        % plateau cross-check, same spirit as pol6's Phase 1, so it does
        % not need pol6's larger per-dim budget -- a single dimension is
        % cheap regardless. 2000 points is ample to confirm the
        % flat-plateau signature and locate the corner to grid resolution
        % as an independent check on the closed form.
        N_GRID = 2000;
    elseif strcmp(pol_id, 'pol6')
        N_GRID = 636 * ones(1, nfree);
    else
        switch pol_id
            case 'pol4', TARGET_TOTAL_GRID = 2197000;   % 130/dim
            case 'pol3', TARGET_TOTAL_GRID = 125000;    % ~354/dim (nfree=2),
                                                         % a conservative
                                                         % default, no real
                                                         % wall-clock timing
                                                         % yet.
            otherwise,   TARGET_TOTAL_GRID = 125000;    % unrecognized pol_id fallback,
                                                         % should not be reachable given
                                                         % the pol_id validation earlier
                                                         % in this file.
        end
        n_per_dim = max(4, round(TARGET_TOTAL_GRID ^ (1 / nfree)));
        N_GRID = n_per_dim * ones(1, nfree);
    end
elseif isscalar(N_GRID)
    N_GRID = N_GRID * ones(1, nfree);
end
if isempty(N_REFINE)
    % Phase 2 (simulated annealing) candidate counts per pol. Phase 2 is
    % parallelized (parfor); wall-clock impact scales with available cores.
    switch pol_id
        case 'pol2', N_REFINE = 60;   % 1-D but genuine interior tradeoff
                                       % (no flat plateau, unlike pol1) --
                                       % more candidates than pol1's 20,
                                       % less than the multi-D pols since
                                       % it's still only 1 dimension.
        case 'pol1', N_REFINE = 20;   % 1-D, closed form already known --
                                       % SA here is a cross-check, not the
                                       % source of truth. Small candidate
                                       % set is sufficient.
        case 'pol6', N_REFINE = 40;
        case 'pol4', N_REFINE = 160;
        case 'pol3', N_REFINE = 100;   % a conservative default, no real
                                        % timing yet to tune against.
        otherwise,   N_REFINE = 100;   % unrecognized pol_id fallback
    end
end


sa_opts = optimoptions('simulannealbnd', ...
    'Display',               'off', ...
    'MaxIterations',          1500 + 500*(nfree-2), ...
    'MaxFunctionEvaluations', 3000 + 1000*(nfree-2), ...
    'InitialTemperature',     0.05, ...
    'ReannealInterval',       200, ...
    'TemperatureFcn',        'temperatureexp', ...
    'AnnealingFcn',          'annealingfast', ...
    'ObjectiveLimit',        -Inf);

% =========================================================================
% BENCHMARK  (tau_R = 0, everything else at its calibrated baseline)
% =========================================================================
x_bench = zeros(1, nfree);
for i = 1:nfree
    switch cfg.free{i}
        case 'tau_R', x_bench(i) = 0.0;
        case 'tau_K', x_bench(i) = fp.tau_K;
        case 'tau_L', x_bench(i) = fp.tau_L;
        case 'w_min', x_bench(i) = w_min_calib;
    end
end
welf_bench = eval_policy_pol(x_bench, cfg, cp);
if ~welf_bench.converged
    error(['welfare_search_pol: BGP failed at the benchmark point. ' ...
           'Check %s is consistent with the current solvers.'], cal_file);
end
SWF_bench = welf_bench.SWF;
fprintf('\nBenchmark SWF = %.8f   regime=%s   u_a1=%.4f%%  u_a2=%.4f%%\n\n', ...
        SWF_bench, welf_bench.regime, 100*welf_bench.u_a1, 100*welf_bench.u_a2);

% =========================================================================
% PHASE 1: N-D GRID SCAN (parallel)
% =========================================================================
fprintf('=== Phase 1: %d-D grid scan (%s = %d points) ===\n', ...
        nfree, strjoin(cellfun(@num2str, num2cell(N_GRID), 'UniformOutput', false), ' x '), ...
        prod(N_GRID));
t_start = tic;

vecs = cell(1, nfree);
for i = 1:nfree
    vecs{i} = linspace(lb(i), ub(i), N_GRID(i))';
end
grids = cell(1, nfree);
[grids{:}] = ndgrid(vecs{:});
N_TOTAL = numel(grids{1});
X_flat = zeros(N_TOTAL, nfree);
for i = 1:nfree
    X_flat(:, i) = grids{i}(:);
end

SWF_flat = -inf(N_TOTAL, 1);
reg_flat = repmat({'?'}, N_TOTAL, 1);

parfor i = 1:N_TOTAL
    w = eval_policy_pol(X_flat(i, :), cfg, cp);   %#ok<PFBNS>
    reg_flat{i} = w.regime;
    if w.converged
        SWF_flat(i) = w.SWF;
    end
end

elapsed1 = toc(t_start);
n_fail1  = sum(~isfinite(SWF_flat));
[SWF_desc, sort_idx] = sort(SWF_flat, 'descend');

fprintf('Phase 1 complete (%.1f s). Failures: %d / %d\n', elapsed1, n_fail1, N_TOTAL);
regime_names = {'competitive', 'mw_binding', 'joint_binding'};
for r = 1:numel(regime_names)
    fprintf('  %-14s : %d\n', regime_names{r}, sum(strcmp(reg_flat, regime_names{r})));
end
fprintf('Best grid SWF = %.8f  at  x = [%s]  (regime: %s)\n\n', ...
        SWF_desc(1), num2str(X_flat(sort_idx(1), :), '%.4f  '), reg_flat{sort_idx(1)});

% =========================================================================
% PHASE 2: SA REFINEMENT (parallel)
% =========================================================================
fprintf('=== Phase 2: %d-D SA refinement (%d runs, parallel) ===\n', nfree, N_REFINE);
t_start = tic;

n_cands = min(N_REFINE, sum(isfinite(SWF_desc)));
X0_all   = X_flat(sort_idx(1:n_cands), :);
SWF0_all = SWF_desc(1:n_cands);

X_ref   = nan(n_cands, nfree);
SWF_ref = -inf(n_cands, 1);

% Simulated annealing draws from the worker's own stream, which is not
% tied to the client's, so the refinement was the one step a rerun could
% not reproduce. Pinning a substream to the candidate index fixes that
% whatever the pool size, since substream j then belongs to candidate j
% rather than to whichever worker happened to take it.
SA_SEED = 42;

parfor j = 1:n_cands
    sa_stream = RandStream('mrg32k3a', 'Seed', SA_SEED);
    sa_stream.Substream = j;
    RandStream.setGlobalStream(sa_stream);
    x0_j = X0_all(j, :);   %#ok<PFBNS>
    try
        [x_opt, neg_opt] = simulannealbnd( ...
            @(x) neg_swf_obj_pol(x, cfg, cp), ...   %#ok<PFBNS>
            x0_j, lb, ub, sa_opts);   %#ok<PFBNS>
        X_ref(j, :) = x_opt;
        SWF_ref(j)  = -neg_opt;
    catch
        X_ref(j, :) = x0_j;
        SWF_ref(j)  = SWF0_all(j);   %#ok<PFBNS>
    end
end

elapsed2 = toc(t_start);
fprintf('Phase 2 complete (%.1f s).\n\n', elapsed2);

[SWF_star, best_j] = max(SWF_ref);
x_star = X_ref(best_j, :);

fprintf('=== GLOBAL OPTIMUM (SA cross-check -- see locate_corner_pol for the exact one) ===\n');
results = struct();
for i = 1:nfree
    fn = [cfg.free{i} '_star'];
    results.(fn) = x_star(i);
    fprintf('  %-6s* = %.6f\n', cfg.free{i}, x_star(i));
end
results.SWF_star = SWF_star;
results.SWF_bench = SWF_bench;
fprintf('  SWF*       = %.8f\n', SWF_star);
fprintf('  SWF(bench) = %.8f\n', SWF_bench);
fprintf('  Delta SWF  = %+.8f\n\n', SWF_star - SWF_bench);

results.SPEC = SPEC;
results.cp   = cp;
results.lb   = lb;
results.ub   = ub;
results.X_flat   = X_flat;     % needed by plot_policy_grid_pol.m to
results.SWF_flat = SWF_flat;   % reconstruct the N-D grid -- computed above
results.reg_flat = reg_flat;   % but was not being persisted; fixed here.
% Baseline value of ALL FOUR instruments, regardless of which are free for
% this pol. Needed by plot_policy_grid_pol.m: a FREE instrument (e.g.
% tau_K for pol4) is never stored in cp, so cp alone cannot answer "what
% is tau_K's calibrated baseline" once the search finishes -- only the
% search itself (here) has fp in scope to answer that.
results.baseline_tau_K = fp.tau_K;
results.baseline_tau_L = fp.tau_L;
results.baseline_w_min = w_min_calib;
save(res_file, '-struct', 'results');
fprintf('Saved %s\n\n', res_file);

end  % ---- main function --------------------------------------------------


% =========================================================================
% LOCAL HELPERS
% =========================================================================
function neg_swf = neg_swf_obj_pol(x, cfg, cp)
    PENALTY = 1e10;
    w = eval_policy_pol(x, cfg, cp);
    if w.converged
        neg_swf = -w.SWF;
    else
        neg_swf = PENALTY;
    end
end

function s = tern(cond, s_true, s_false)
    if cond, s = s_true; else, s = s_false; end
end
