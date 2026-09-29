function verify_manuscript_labels()
% VERIFY_MANUSCRIPT_LABELS  Checks the policy registry still says
% what it should. Worth re-running after any edit near pol_config or
% manuscript_policy_map.
%
%   verify_manuscript_labels()
%
% cfg.tag feeds pol_filenames and so the names of result files already
% computed and checked; a stray edit there would misfile results or hide
% them behind a name nothing looks for. Three checks, stopping at the first
% failure:
%
%   1. cfg.tag and cfg.free, rebuilt through pol_config, match the strings
%      recorded in EXPECTED_TAG and EXPECTED_FREE below.
%   2. manuscript_policy_map is one-to-one onto the Policy numbers the
%      registered identifiers carry, with none mapped twice.
%   3. The label pol_config attaches agrees with what
%      manuscript_policy_map returns.
%
% Runs under plain Octave, needs no toolboxes, touches no .mat files.
POL_IDS = {'pol1','pol2','pol3','pol4','pol5','pol6'};

% The policy identifier and the manuscript number coincide for all five.
EXPECTED_TAG = struct( ...
    'pol1', 'pol1_tauR', ...
    'pol2', 'pol2_tauK', ...
    'pol3', 'pol3_tauK_tauL', ...
    'pol4', 'pol4_tauR_tauK_tauL', ...
    'pol5', 'pol5_wmin', ...
    'pol6', 'pol6_tauR_wmin');

EXPECTED_FREE = struct( ...
    'pol1', {{'tau_R'}}, ...
    'pol2', {{'tau_K'}}, ...
    'pol3', {{'tau_K','tau_L'}}, ...
    'pol4', {{'tau_R','tau_K','tau_L'}}, ...
    'pol5', {{'w_min'}}, ...
    'pol6', {{'tau_R','w_min'}});

fprintf('--- verify_manuscript_labels: checking %d pols ---\n', numel(POL_IDS));

% Check 1: the identifier, tag and free set match what was recorded
for i = 1:numel(POL_IDS)
    pid = POL_IDS{i};
    cfg = pol_config(pid);

    if ~strcmp(cfg.id, pid)
        error('verify_manuscript_labels: %s: cfg.id=''%s'' (expected ''%s'').', ...
              pid, cfg.id, pid);
    end
    if ~strcmp(cfg.tag, EXPECTED_TAG.(pid))
        error(['verify_manuscript_labels: %s: cfg.tag=''%s'' but expected ' ...
               '''%s'' -- this would misfile .mat results.'], ...
               pid, cfg.tag, EXPECTED_TAG.(pid));
    end
    exp_free = EXPECTED_FREE.(pid);
    if numel(cfg.free) ~= numel(exp_free) || ~isequal(cfg.free, exp_free)
        error('verify_manuscript_labels: %s: cfg.free changed from expected set.', pid);
    end
    if ~isfield(cfg, 'manuscript_label') || isempty(cfg.manuscript_label)
        error('verify_manuscript_labels: %s: cfg.manuscript_label missing.', pid);
    end
    if ~strcmp(cfg.manuscript_label, sprintf('Policy %s', pid(4:end)))
        error(['verify_manuscript_labels: %s: cfg.manuscript_label=''%s'' is ' ...
               'not the identity label expected post-renumbering.'], ...
               pid, cfg.manuscript_label);
    end
    fprintf('  [ok] %-5s  tag=%-22s  free={%s}%-4s  -> %s\n', ...
            pid, cfg.tag, strjoin(cfg.free, ','), '', cfg.manuscript_label);
end

% Check 2: the map runs one-to-one onto the registered Policy numbers
labels = cell(1, numel(POL_IDS));
for i = 1:numel(POL_IDS)
    labels{i} = manuscript_policy_map(POL_IDS{i});
end
expected_labels = cellfun(@(p) sprintf('Policy %s', p(4:end)), POL_IDS, 'UniformOutput', false);
if numel(unique(labels)) ~= numel(POL_IDS)
    error('verify_manuscript_labels: manuscript_policy_map.m is NOT injective -- duplicate label found.');
end
if ~isempty(setxor(labels, expected_labels))
    error(['verify_manuscript_labels: manuscript_policy_map.m does not cover ' ...
           'the registered policies exactly once.']);
end

% Check 3: the label pol_config attaches agrees with the map
for i = 1:numel(POL_IDS)
    pid = POL_IDS{i};
    cfg = pol_config(pid);
    want = manuscript_policy_map(pid);
    if ~strcmp(cfg.manuscript_label, want)
        error(['verify_manuscript_labels: %s: pol_config.m''s ' ...
               'cfg.manuscript_label=''%s'' disagrees with ' ...
               'manuscript_policy_map(''%s'')=''%s''.'], ...
               pid, cfg.manuscript_label, pid, want);
    end
end

fprintf('\nALL CHECKS PASSED -- every registered pol coincides with the\n');
fprintf('manuscript''s own Policy numbering; cfg.tag/.free match the\n');
fprintf('recorded set for all %d pols; manuscript_policy_map.m is a\n', numel(POL_IDS));
fprintf('verified identity bijection.\n');

end
