function sc = compute_scalars(params)
% COMPUTE_SCALARS  Every derived scalar the BGP solver needs.
%
%   sc = compute_scalars(params)
%

alpha   = params.alpha;
psi     = params.psi;
theta   = params.theta;
eta_bar = params.eta_bar;
lam_bar = params.lam_bar;
r       = params.r;
delta   = params.delta;
tau_R   = params.tau_R;
s_n     = params.s_n;
s_a1    = params.s_a1;
s_a2    = params.s_a2;

% -------------------------------------------------------------------------
% Check the parameters make sense
% -------------------------------------------------------------------------
assert(alpha > 0 && alpha < 1,        'alpha must be in (0,1)');
assert(psi   > 0 && psi   < 1,        'psi must be in (0,1)');
assert(alpha + psi < 1,               'alpha + psi must be < 1');
assert(theta < 0,                     'theta must be < 0 (tasks are complements)');
assert(eta_bar > 1,                   'eta_bar must be > 1');
assert(lam_bar > 0,                   'lam_bar must be > 0');
assert(r > 0,                         'r must be positive');
assert(delta >= 0,                    'delta must be non-negative');
assert(tau_R > -1,                    'tau_R must exceed -1 (so that p = (1+tau_R)*(r+delta) > 0)');
assert(s_n  > 0,                      's_n must be positive');
assert(s_a1 > 0,                      's_a1 must be positive');
assert(s_a2 > 0,                      's_a2 must be positive');

% -------------------------------------------------------------------------
% The derived scalars
% -------------------------------------------------------------------------

% The rental price of robots
p = (1 + tau_R) * (r + delta);

% Joint density of (eta, lambda) ~ Uniform[1,eta_bar] x Uniform[0,lam_bar]
f = 1 / (lam_bar * (eta_bar - 1));

gam = theta / (1 - theta);
g1  = gam + 1;   % between 0 and 1; positive, so the power integrals stay finite
g2  = gam + 2;   % between 1 and 2, and positive for the same reason

% Production function exponents
exp_y     =  psi / (1 - alpha);               % y = A * z^exp_y
exp_z     = (1 - alpha - psi) / (1 - alpha);  % z-exponent in (E1); > 0
exp_Omega = (1 - theta) / theta;              % Omega-exponent in (E1); < 0

A = (alpha / (r + delta))^(alpha / (1 - alpha)) ...
    * s_n^((1 - alpha - psi) / (1 - alpha));

% -------------------------------------------------------------------------
% Assemble the output, carrying the inputs along for convenience
% -------------------------------------------------------------------------
sc.alpha     = alpha;
sc.psi       = psi;
sc.theta     = theta;
sc.eta_bar   = eta_bar;
sc.lam_bar   = lam_bar;
sc.r         = r;
sc.delta     = delta;
sc.tau_R     = tau_R;
sc.s_n       = s_n;
sc.s_a1      = s_a1;
sc.s_a2      = s_a2;

sc.p         = p;
sc.f         = f;
sc.gamma     = gam;
sc.g1        = g1;
sc.g2        = g2;
sc.A         = A;
sc.exp_y     = exp_y;
sc.exp_z     = exp_z;
sc.exp_Omega = exp_Omega;

end
