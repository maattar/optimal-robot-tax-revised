function cfg = pol_config(pol_id)
% POL_CONFIG  policy instruments 
%
%   pol1  tau_R                   pol4  tau_R, tau_K, tau_L
%   pol2  tau_K                   pol5  w_min
%   pol3  tau_K, tau_L            pol6  tau_R, w_min
%
% 
switch pol_id
    case 'pol1'
        cfg.id          = 'pol1';
        cfg.free        = {'tau_R'};
        cfg.needs_euler = false;
        cfg.tag         = 'pol1_tauR';
        cfg.built       = true;

    case 'pol2'
        cfg.id          = 'pol2';
        cfg.free        = {'tau_K'};
        cfg.needs_euler = true;
        cfg.tag         = 'pol2_tauK';
        cfg.built       = true;

    case 'pol3'
        cfg.id          = 'pol3';
        cfg.free        = {'tau_K', 'tau_L'};
        cfg.needs_euler = true;
        cfg.tag         = 'pol3_tauK_tauL';
        cfg.built       = true;

    case 'pol4'
        cfg.id          = 'pol4';
        cfg.free        = {'tau_R', 'tau_K', 'tau_L'};
        cfg.needs_euler = true;
        cfg.tag         = 'pol4_tauR_tauK_tauL';
        cfg.built       = true;

    case 'pol5'
        cfg.id          = 'pol5';
        cfg.free        = {'w_min'};
        cfg.needs_euler = false;
        cfg.tag         = 'pol5_wmin';
        cfg.built       = true;

    case 'pol6'
        cfg.id          = 'pol6';
        cfg.free        = {'tau_R', 'w_min'};
        cfg.needs_euler = false;
        cfg.tag         = 'pol6_tauR_wmin';
        cfg.built       = true;

    otherwise
        error(['pol_config: unrecognized pol_id ''%s''. ' ...
               'Valid: pol1, pol2, pol3, pol4, pol5, pol6.'], pol_id);
end


cfg.manuscript_label = manuscript_policy_map(cfg.id);

end
