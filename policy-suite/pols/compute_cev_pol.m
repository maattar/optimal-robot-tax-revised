function cev = compute_cev_pol(rep0, rep1, m_sq, m_opt)
% COMPUTE_CEV_POL  Consumption equivalent variations between the status quo
% path and an optimal one
GROUPS = {'U_a1'; 'U_a2'; 'E_a1'; 'E_a2'; 'E_n'; 'Q3'; 'Q4'; 'Q5'};


POP_EPS = 1e-9;

if nargin < 3, m_sq  = []; end
if nargin < 4, m_opt = []; end

cev = struct();
cev.groups = GROUPS;
cev.sigma  = local_sigma(rep0, rep1);

if ~local_report_ok(rep0) || ~local_report_ok(rep1)
    cev = local_fill_nan(cev);
    return
end

c_sq  = local_cvec(rep0);
c_opt = local_cvec(rep1);

if isempty(m_sq)
    m_sq = local_mvec(rep0);
else
    m_sq = m_sq(:);
end
if isempty(m_opt)
    m_opt = local_mvec(rep1);
else
    m_opt = m_opt(:);
end

n = numel(GROUPS);
if numel(m_sq) ~= n || numel(m_opt) ~= n
    error(['compute_cev_pol: population share vectors must have %d entries ' ...
           '(got %d and %d), in the order %s.'], ...
          n, numel(m_sq), numel(m_opt), strjoin(GROUPS', ', '));
end

cev.c_sq  = c_sq;   cev.c_opt  = c_opt;
cev.m_sq  = m_sq;   cev.m_opt  = m_opt;

sigma = cev.sigma;

% By group
CEV_g = c_opt ./ c_sq - 1;
vanishing = (m_sq <= POP_EPS) | (m_opt <= POP_EPS) | ~isfinite(c_sq) | (c_sq == 0);
CEV_g(vanishing) = NaN;

cev.CEV_g     = CEV_g;
cev.CEV_g_pct = 100 * CEV_g;

% The reallocatd workers
REALLOC = {'a1'; 'a2'};
iU = [1; 2];
iE = [3; 4];
nr = numel(REALLOC);

realloc_mass  = zeros(nr, 1);
realloc_c_sq  = NaN(nr, 1);
realloc_c_opt = NaN(nr, 1);

for i = 1:nr
    realloc_mass(i) = m_sq(iU(i)) - m_opt(iU(i));
    if realloc_mass(i) > POP_EPS
        realloc_c_sq(i)  = c_sq(iU(i));
        realloc_c_opt(i) = c_opt(iE(i));
    elseif realloc_mass(i) < -POP_EPS
        realloc_c_sq(i)  = c_sq(iE(i));
        realloc_c_opt(i) = c_opt(iU(i));
    end
end

realloc_CEV = realloc_c_opt ./ realloc_c_sq - 1;
realloc_CEV(~isfinite(realloc_c_sq) | realloc_c_sq <= 0) = NaN;

cev.realloc_types   = REALLOC;
cev.realloc_mass    = realloc_mass;
cev.realloc_c_sq    = realloc_c_sq;
cev.realloc_c_opt   = realloc_c_opt;
cev.realloc_CEV     = realloc_CEV;
cev.realloc_CEV_pct = 100 * realloc_CEV;

exante_CEV    = NaN(nr, 1);
exante_Ec_sq  = NaN(nr, 1);
exante_Ec_opt = NaN(nr, 1);

for i = 1:nr
    idx = [iU(i); iE(i)];
    [S0, Ec0] = local_lottery(m_sq,  c_sq,  idx, sigma);
    [S1, Ec1] = local_lottery(m_opt, c_opt, idx, sigma);
    exante_Ec_sq(i)  = Ec0;
    exante_Ec_opt(i) = Ec1;
    if ~isfinite(S0) || ~isfinite(S1)
        continue
    end
    if abs(1 - sigma) < 1e-12
        exante_CEV(i) = exp(S1 - S0) - 1;
    elseif S0 > 0 && S1 > 0
        exante_CEV(i) = (S1 / S0) ^ (1 / (1 - sigma)) - 1;
    end
end

cev.exante_CEV     = exante_CEV;
cev.exante_CEV_pct = 100 * exante_CEV;
cev.exante_Ec_sq   = exante_Ec_sq;
cev.exante_Ec_opt  = exante_Ec_opt;

% For the economy as a whole
a_sq  = m_sq  > 0;
a_opt = m_opt > 0;

cev.active_sq  = a_sq;
cev.active_opt = a_opt;
cev.mass_sq    = sum(m_sq(a_sq));
cev.mass_opt   = sum(m_opt(a_opt));

bad = any(~isfinite(c_sq(a_sq)))  || any(c_sq(a_sq)  <= 0) || ...
      any(~isfinite(c_opt(a_opt))) || any(c_opt(a_opt) <= 0) || ...
      ~isfinite(sigma);

if bad
    cev.S_sq = NaN; cev.S_opt = NaN;
    cev.CEV_agg = NaN; cev.CEV_agg_pct = NaN;
    cev.CEV_agg_from_SWF = NaN; cev.identity_gap = NaN;
    cev.converged = false;
    return
end

if abs(1 - sigma) < 1e-12
    
    S_sq  = sum(m_sq(a_sq)   .* log(c_sq(a_sq)));
    S_opt = sum(m_opt(a_opt) .* log(c_opt(a_opt)));
    CEV_agg = exp(S_opt - S_sq) - 1;
    CEV_agg_from_SWF = NaN;   % the identity degenerates at this point
else
    S_sq  = sum(m_sq(a_sq)   .* c_sq(a_sq)   .^ (1 - sigma));
    S_opt = sum(m_opt(a_opt) .* c_opt(a_opt) .^ (1 - sigma));
    if S_sq <= 0 || S_opt <= 0
        CEV_agg = NaN;
    else
        CEV_agg = (S_opt / S_sq) ^ (1 / (1 - sigma)) - 1;
    end
    CEV_agg_from_SWF = local_cev_from_swf(rep0.SWF, rep1.SWF, sigma);
end

cev.S_sq             = S_sq;
cev.S_opt            = S_opt;
cev.CEV_agg          = CEV_agg;
cev.CEV_agg_pct      = 100 * CEV_agg;
cev.CEV_agg_from_SWF = CEV_agg_from_SWF;
if isfinite(CEV_agg) && isfinite(CEV_agg_from_SWF)
    cev.identity_gap = abs(CEV_agg_from_SWF - CEV_agg);
else
    cev.identity_gap = NaN;
end
cev.converged = isfinite(CEV_agg);

end  


function sigma = local_sigma(rep0, rep1)
  
    if ~isfield(rep0, 'sigma') || ~isfield(rep1, 'sigma')
        error(['compute_cev_pol: report struct has no .sigma field. ' ...
               'CEV^agg cannot be computed without the CRRA parameter. ' ...
               'This usually means the reports were produced by an older ' ...
               'compute_policy_report_pol.m -- regenerate them.']);
    end
    s0 = rep0.sigma;  s1 = rep1.sigma;
    if isfinite(s0) && isfinite(s1) && abs(s0 - s1) > 1e-12
        error(['compute_cev_pol: the two reports carry different sigma ' ...
               '(%.12g vs %.12g). They are not comparable -- a CEV between ' ...
               'them would be meaningless.'], s0, s1);
    end
    if isfinite(s0)
        sigma = s0;
    else
        sigma = s1;
    end
end

function ok = local_report_ok(rep)
    ok = isfield(rep, 'converged') && rep.converged;
end

function c = local_cvec(rep)
    c = [rep.c_Ua1; rep.c_Ua2; rep.c_Ea1; rep.c_Ea2; ...
         rep.c_En;  rep.c_Q3;  rep.c_Q4;  rep.c_Q5];
end

function m = local_mvec(rep)
  

    n_En = 0.4 - (rep.u_a1 + rep.L_a1) - (rep.u_a2 + rep.L_a2);
    m = [rep.u_a1; rep.u_a2; rep.L_a1; rep.L_a2; n_En; 0.2; 0.2; 0.2];
end

function cev_val = local_cev_from_swf(SWF0, SWF1, sigma)
    
    S0 = 1 + (1 - sigma) * SWF0;
    S1 = 1 + (1 - sigma) * SWF1;
    if ~isfinite(S0) || ~isfinite(S1) || S0 <= 0 || S1 <= 0
        cev_val = NaN;
        return
    end
    cev_val = (S1 / S0) ^ (1 / (1 - sigma)) - 1;
end

function cev = local_fill_nan(cev)
    n = numel(cev.groups);
    cev.c_sq = NaN(n,1);  cev.c_opt = NaN(n,1);
    cev.m_sq = NaN(n,1);  cev.m_opt = NaN(n,1);
    cev.active_sq = false(n,1);  cev.active_opt = false(n,1);
    cev.mass_sq = NaN;  cev.mass_opt = NaN;
    cev.CEV_g = NaN(n,1);  cev.CEV_g_pct = NaN(n,1);
    cev.realloc_types = {'a1'; 'a2'};
    cev.realloc_mass = NaN(2,1);
    cev.realloc_c_sq = NaN(2,1);  cev.realloc_c_opt = NaN(2,1);
    cev.realloc_CEV = NaN(2,1);  cev.realloc_CEV_pct = NaN(2,1);
    cev.exante_CEV = NaN(2,1);  cev.exante_CEV_pct = NaN(2,1);
    cev.exante_Ec_sq = NaN(2,1);  cev.exante_Ec_opt = NaN(2,1);
    cev.S_sq = NaN;  cev.S_opt = NaN;
    cev.CEV_agg = NaN;  cev.CEV_agg_pct = NaN;
    cev.CEV_agg_from_SWF = NaN;  cev.identity_gap = NaN;
    cev.converged = false;
end

function [S, Ec] = local_lottery(m, c, idx, sigma)
    
    mm = m(idx);  cc = c(idx);
    act = mm > 0;
    s_tot = sum(mm);
    if s_tot <= 0 || any(~isfinite(cc(act))) || any(cc(act) <= 0)
        S = NaN;  Ec = NaN;
        return
    end
    if abs(1 - sigma) < 1e-12
        S = sum(mm(act) .* log(cc(act))) / s_tot;
    else
        S = sum(mm(act) .* cc(act) .^ (1 - sigma)) / s_tot;
    end
    Ec = sum(mm(act) .* cc(act)) / s_tot;
end
