function result = check_wmin_identification_pol(pol_id, SPEC, x_override)
% CHECK_WMIN_IDENTIFICATION_POL  For POL5 and POL6, check whether the 
% minwage is a real optimum, or one arbitrary point on a flat surface
if nargin < 2 || isempty(SPEC), SPEC = 'baseline'; end
if nargin < 3,                  x_override = [];   end

if strcmp(pol_id, 'pol4')
    error(['check_wmin_identification_pol: pol4 holds w_min FIXED at ' ...
           'the calibrated baseline for its whole search (see ' ...
           'pol_config.m) -- there is no w_min* to interrogate. pol5 and ' ...
           'pol6 are the policies that free the floor.']);
end
if ~any(strcmp(pol_id, {'pol5', 'pol6'}))
    error(['check_wmin_identification_pol: unrecognized/unsupported ' ...
           'pol_id ''%s''. Valid: pol5, pol6.'], pol_id);
end

cfg = pol_config(pol_id);
if ~cfg.built
    error('check_wmin_identification_pol: %s is not marked built in pol_config.m.', pol_id);
end
if ~any(strcmp('w_min', cfg.free))
    error(['check_wmin_identification_pol: %s does not free w_min ' ...
           '(cfg.free = %s) -- nothing to check.'], pol_id, strjoin(cfg.free, ', '));
end

[cal_file, res_file] = pol_filenames(cfg, SPEC);
out_file = wmin_id_filename(cfg, SPEC);

% =========================================================================
% LOAD CALIB
if ~exist(cal_file, 'file')
    error('check_wmin_identification_pol: %s not found.', cal_file);
end
S = load(cal_file);
if ~isfield(S, 'zeta_star') || ~isfield(S, 'fp')
    error('check_wmin_identification_pol: %s must contain both ''zeta_star'' and ''fp''.', cal_file);
end
if isfield(S, 'SPEC') && ~isempty(SPEC) && ~strcmp(S.SPEC, SPEC)
    error(['check_wmin_identification_pol: SPEC MISMATCH. Requested ' ...
           '''%s'' but %s identifies itself as ''%s''.'], SPEC, cal_file, S.SPEC);
end
zeta_star = S.zeta_star;  fp = S.fp;
w_min_calib = zeta_star(5);

cp.alpha   = zeta_star(1);
cp.psi     = zeta_star(2);
cp.lam_bar = zeta_star(3);
cp.eta_bar = zeta_star(4);
cp.delta = fp.delta;  cp.theta = fp.theta;  cp.nu   = fp.nu;
cp.sigma = fp.sigma;  cp.g_A   = fp.g_A;    cp.g_N  = fp.g_N;
cp.rho   = fp.rho;    cp.s_a1  = fp.s_a1;   cp.s_a2 = fp.s_a2;
cp.s_n   = 1 - fp.s_a1 - fp.s_a2;
cp.verbose = false;

fprintf('\n===========================================================\n');
fprintf('  CHECK_WMIN_IDENTIFICATION_POL(''%s'')  SPEC=''%s''\n', pol_id, SPEC);
fprintf('===========================================================\n');
fprintf('  free = %s\n', strjoin(cfg.free, ', '));

if cfg.needs_euler
    [Phi, ~] = phi_from_calibration(fp, false);
    cp.Phi   = Phi;
    r_base_recon = euler_closure(Phi, fp.tau_K);
    fprintf('  Euler round-trip r(tau_K_base)==fp.r gap = %.3e\n', abs(r_base_recon - fp.r));
else
    cp.r = fp.r;
end

if ~any(strcmp('tau_K', cfg.free)), cp.tau_K = fp.tau_K;    end
if ~any(strcmp('tau_L', cfg.free)), cp.tau_L = fp.tau_L;    end
if ~any(strcmp('tau_R', cfg.free)), cp.tau_R = 0;           end

nfree     = numel(cfg.free);
wmin_idx  = find(strcmp(cfg.free, 'w_min'));

if ~isempty(x_override)
    if numel(x_override) ~= nfree
        error(['check_wmin_identification_pol: x_override has %d ' ...
               'element(s) but %s expects %d (%s).'], numel(x_override), ...
              pol_id, nfree, strjoin(cfg.free, ', '));
    end
    x_star = x_override(:)';
    x_src  = 'explicit x_override';
else
    if ~exist(res_file, 'file')
        error(['check_wmin_identification_pol: %s not found, and no ' ...
               'x_override given. Run welfare_search_pol(''%s'', ''%s'') ' ...
               'first, or pass the point to interrogate explicitly.'], ...
              res_file, pol_id, SPEC);
    end
    R = load(res_file);
    x_star = zeros(1, nfree);
    missing = {};
    for i = 1:nfree
        fn = [cfg.free{i} '_star'];
        if isfield(R, fn)
            x_star(i) = R.(fn);
        else
            missing{end+1} = fn; %#ok<AGROW>
        end
    end
    if ~isempty(missing)
        error(['check_wmin_identification_pol: %s is missing field(s): ' ...
               '%s.'], res_file, strjoin(missing, ', '));
    end
    x_src = sprintf('welfare_result (%s)', res_file);
end

% PLATEAU_EPS 
PLATEAU_EPS = 1e-6;     
tauR_idx = find(strcmp(cfg.free, 'tau_R'));
if ~isempty(tauR_idx)
    x_star(tauR_idx) = x_star(tauR_idx) + PLATEAU_EPS;
    fprintf('  [PLATEAU_EPS guard: tau_R nudged by +%.1e off the corner, project convention]\n', ...
            PLATEAU_EPS);
end

w_min_star = x_star(wmin_idx);
fprintf('\n  x_star = [%s]   (source: %s)\n', num2str(x_star, '%.6f  '), x_src);
fprintf('  w_min* = %.9f\n\n', w_min_star);





% THRESHOLD
x_thresh = x_star;  x_thresh(wmin_idx) = 0;
welf0 = eval_policy_pol(x_thresh, cfg, cp);
if ~welf0.converged || ~strcmp(welf0.regime, 'competitive')
    error(['check_wmin_identification_pol: evaluating at w_min=0 (all ' ...
           'other instruments held at x_star) did not return a converged ' ...
           '''competitive'' regime (got regime=''%s'', converged=%d) -- ' ...
           'the threshold wage cannot be identified. This should be ' ...
           'impossible (w_min=0 is <= any positive wage); investigate ' ...
           'before trusting anything else in this report.'], ...
          welf0.regime, welf0.converged);
end
threshold = welf0.w_a1;

welf_star = eval_policy_pol(x_star, cfg, cp);
fprintf('--- Threshold and regime at x_star ---\n');
fprintf('  threshold (competitive w_a1 at w_min=0) = %.12f\n', threshold);
fprintf('  regime AT x_star (w_min=%.6f)            = %s   (converged=%d)\n\n', ...
        w_min_star, tern2(welf_star.converged, welf_star.regime, 'failed'), welf_star.converged);


% projection
TOL_FLAT     = 1e-9;
N_SUBGRID    = 41;

sub_grid = linspace(0, threshold * (1 - 1e-6), N_SUBGRID);
sub_SWF  = -inf(1, N_SUBGRID);
sub_reg  = repmat({'?'}, 1, N_SUBGRID);
ws = warning('off', 'all');   


for i = 1:N_SUBGRID
    x = x_star; x(wmin_idx) = sub_grid(i);
    w = eval_policy_pol(x, cfg, cp);
    sub_reg{i} = w.regime;
    if w.converged, sub_SWF(i) = w.SWF; end
end
warning(ws);
flat_gap = max(abs(sub_SWF - sub_SWF(1)));
n_not_competitive = sum(~strcmp(sub_reg, 'competitive'));

fprintf('--- Sub-threshold flatness sweep: w_min in [0, threshold), %d points ---\n', N_SUBGRID);
fprintf('  max|SWF-SWF(0)| = %.3e   %s\n', flat_gap, tern2(flat_gap < TOL_FLAT, ...
        'FLAT to machine precision (as theory predicts)', '**** NOT FLAT -- unexpected'));
if n_not_competitive > 0
    fprintf('  **** %d/%d sub-threshold grid points did NOT return regime=competitive -- unexpected\n', ...
            n_not_competitive, N_SUBGRID);
end
fprintf('\n');

% Local check
below_threshold = (w_min_star < threshold - PLATEAU_EPS);
at_threshold     = (abs(w_min_star - threshold) <= PLATEAU_EPS);

if below_threshold
    local_grid = sub_grid;  local_SWF = sub_SWF;
    local_note = 'w_min* lies strictly below threshold -- covered by the sub-threshold sweep above.';
else
    N_LOCAL = 5;
    if at_threshold
        span = max(1e-4, 0.02 * max(threshold, 1e-6));
        local_grid = w_min_star + span * linspace(-1, 1, N_LOCAL);
        local_grid = max(local_grid, 0);
    else
      span = min(0.25 * (w_min_star - threshold), 0.02 * max(w_min_star, 1e-6));
        span = max(span, 1e-6);
        local_grid = w_min_star + span * linspace(-1, 1, N_LOCAL);
        local_grid = max(local_grid, threshold + PLATEAU_EPS/10);
    end
    local_SWF = -inf(1, N_LOCAL);
    local_reg = repmat({'?'}, 1, N_LOCAL);
    for i = 1:N_LOCAL
        x = x_star; x(wmin_idx) = local_grid(i);
        w = eval_policy_pol(x, cfg, cp);
        local_reg{i} = w.regime;
        if w.converged, local_SWF(i) = w.SWF; end
    end
    local_note = sprintf('regimes at local grid: %s', strjoin(local_reg, ', '));
end

fprintf('--- Local check around w_min* ---\n');
fprintf('  %s\n', local_note);
if ~below_threshold
    fprintf('  w_min grid    : %s\n', num2str(local_grid, '%.9f  '));
    fprintf('  SWF           : %s\n\n', num2str(local_SWF, '%.10f  '));
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
is_indeterminate = false; is_boundary_kink = false;
is_genuine        = false; is_anomaly       = false;

if below_threshold
    if flat_gap < TOL_FLAT && n_not_competitive == 0
        is_indeterminate = true;
        verdict = sprintf(['INDETERMINATE. w_min* = %.6f sits strictly below the ' ...
            'threshold wage %.6f, and SWF is exactly flat (max gap %.3e) across ' ...
            'the ENTIRE [0, threshold) interval. w_min* is NOT uniquely ' ...
            'identified -- any w_min in [0, threshold) is welfare-equivalent; ' ...
            'the reported value is an arbitrary point where the search stopped. ' ...
            'Recommended focal point: w_min* = 0.'], w_min_star, threshold, flat_gap);
    else
        is_anomaly = true;
        verdict = sprintf(['ANOMALY. w_min* = %.6f sits below the threshold wage ' ...
            '%.6f, but the sub-threshold sweep is NOT flat (max gap %.3e, ' ...
            '%d/%d points off-regime). Theory predicts exact flatness here -- ' ...
            'this contradicts it. Do not trust w_min* as reported; investigate ' ...
            'for a bug before using this result.'], w_min_star, threshold, flat_gap, ...
            n_not_competitive, N_SUBGRID);
    end
elseif at_threshold
    % Compare the plateau value (sub-threshold SWF) to the just-above-
    % threshold value to see which side is actually better.
    swf_plateau = sub_SWF(1);
    [swf_local_max, j_max] = max(local_SWF);
    rising = (swf_local_max > swf_plateau + TOL_FLAT);
    is_boundary_kink = true;
    if rising
        verdict = sprintf(['BOUNDARY/KINK, RISING. w_min* = %.6f sits at the ' ...
            'regime threshold (%.6f), and SWF strictly INCREASES moving into ' ...
            'the binding regime (plateau SWF=%.10f vs. best local SWF=%.10f at ' ...
            'w_min=%.9f). This suggests the true optimum is INSIDE the binding ' ...
            'region, not at the boundary -- the search likely under-shot; refine ' ...
            'with a finer grid/SA run above threshold=%.6f before treating %.6f ' ...
            'as final.'], w_min_star, threshold, swf_plateau, swf_local_max, ...
            local_grid(j_max), threshold, w_min_star);
    else
        verdict = sprintf(['BOUNDARY/KINK, FLAT-OR-FALLING. w_min* = %.6f sits at ' ...
            'the regime threshold (%.6f); SWF does not rise moving into the ' ...
            'binding regime (plateau SWF=%.10f vs. best local SWF=%.10f). This ' ...
            'is consistent with the SAME indifference-plateau phenomenon as the ' ...
            'below-threshold case, just landed exactly on its right edge -- ' ...
            'treat as INDETERMINATE with focal point w_min*=0, same as the ' ...
            'strictly-below-threshold verdict.'], w_min_star, threshold, ...
            swf_plateau, swf_local_max);
        is_boundary_kink = false;
        is_indeterminate = true;
    end
else
    % Genuinely above threshold: local max test.
    swf_center = welf_star.SWF;
    is_local_max = all(local_SWF <= swf_center + TOL_FLAT);
    [worst_gap, j_worst] = max(local_SWF - swf_center);
    if is_local_max
        is_genuine = true;
        gaps = swf_center - local_SWF;
        gaps(local_grid == w_min_star) = NaN;
        verdict = sprintf(['GENUINE. w_min* = %.6f sits strictly ABOVE the ' ...
            'threshold wage %.6f (floor binds), and SWF is locally maximized ' ...
            'there: SWF(w_min*)=%.10f is >= every nearby evaluation (worst-case ' ...
            'margin %.3e). This is a real, identified optimum, not a plateau ' ...
            'artifact -- it deserves its own economic account, not folding into ' ...
            'an indeterminacy statement.'], w_min_star, threshold, swf_center, ...
            max(gaps(~isnan(gaps)), [], 'omitnan'));
    else
        is_anomaly = true;
        verdict = sprintf(['ANOMALY. w_min* = %.6f sits above the threshold wage ' ...
            '%.6f, but SWF is HIGHER nearby (at w_min=%.9f, SWF=%.10f > ' ...
            'SWF(w_min*)=%.10f by %.3e). The reported w_min* is NOT a local ' ...
            'optimum -- likely a grid-quantization or SA-refinement artifact; ' ...
            'rerun welfare_search_pol with a finer grid/more SA restarts ' ...
            'near this point before trusting it.'], w_min_star, threshold, ...
            local_grid(j_worst), local_SWF(j_worst), swf_center, worst_gap);
    end
end

fprintf('=== VERDICT ===\n%s\n\n', verdict);

% =========================================================================
% SAVE + RETURN
% =========================================================================
result = struct();
result.pol_id           = pol_id;
result.SPEC                = SPEC;
result.x_star               = x_star;
result.x_src                = x_src;
result.w_min_star           = w_min_star;
result.threshold_wage       = threshold;
result.regime_at_xstar      = welf_star.regime;
result.subthreshold_grid    = sub_grid;
result.subthreshold_SWF     = sub_SWF;
result.flat_gap             = flat_gap;
result.local_grid           = local_grid;
result.local_SWF            = local_SWF;
result.verdict               = verdict;
result.is_indeterminate     = is_indeterminate;
result.is_boundary_kink      = is_boundary_kink;
result.is_genuine            = is_genuine;
result.is_anomaly            = is_anomaly;

save(out_file, '-struct', 'result');
fprintf('Saved %s\n\n', out_file);

end 






function fname = wmin_id_filename(cfg, SPEC)
    if isempty(SPEC)
        fname = sprintf('wmin_identification_%s.mat', cfg.tag);
    else
        fname = sprintf('wmin_identification_%s_%s.mat', cfg.tag, SPEC);
    end
end

function s = tern2(cond, a, b)
    if cond, s = a; else, s = b; end
end
