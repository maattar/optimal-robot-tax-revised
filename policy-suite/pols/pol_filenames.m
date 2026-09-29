function [cal_file, res_file, out_file] = pol_filenames(cfg, SPEC)
% POL_FILENAMES  policy and specification
%   SPEC  one of '', baseline, auto33, equal_shares, acemoglu2020.
VALID_SPECS = {'', 'baseline', 'auto33', 'equal_shares', 'acemoglu2020'};
if ~ischar(SPEC) || ~any(strcmp(SPEC, VALID_SPECS))
    error(['pol_filenames: unrecognised SPEC. Valid options: ''''  ' ...
           '(plain file), ''baseline'', ''auto33'', ' ...
           '''equal_shares'', ''acemoglu2020''.']);
end
if ~isfield(cfg, 'tag') || isempty(cfg.tag)
    error('pol_filenames: cfg.tag is required (from pol_config(id)).');
end

if isempty(SPEC)
    cal_file = 'calibration_result.mat';
    res_file = sprintf('welfare_result_%s.mat', cfg.tag);
    out_file = sprintf('corner_result_%s.mat', cfg.tag);
else
    cal_file = sprintf('calibration_result_%s.mat', SPEC);
    res_file = sprintf('welfare_result_%s_%s.mat', cfg.tag, SPEC);
    out_file = sprintf('corner_result_%s_%s.mat', cfg.tag, SPEC);
end

end
