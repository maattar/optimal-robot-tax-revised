function sol = solve_joint_binding_bgp(params, w_min)
% SOLVE_JOINT_BINDING_BGP  
%
%   sol = solve_joint_binding_bgp(params, w_min)
%

sc    = compute_scalars(params);
psi_A = sc.psi * sc.A;
gam   = sc.gamma;


%  Step 1: the regime

hat_lam_fixed = sc.p / w_min;
auto_on       = (hat_lam_fixed < sc.lam_bar);

if ~auto_on
    sol_regime = 'no_auto';
else
    % The split uses the type-2 wage directly, which is the floor.
    if sc.lam_bar < (sc.p / w_min) * sc.eta_bar
        sol_regime = 'case1';
    else
        sol_regime = 'case2';
    end
end

if params.verbose
    fprintf('\n=== Joint-Binding BGP Solver (w_min = %.6f) ===\n', w_min);
    fprintf('Regime (exogenous, no root-find): %s\n', sol_regime);
    fprintf('  hat_lambda = %.6f   lam_bar = %.6f\n', hat_lam_fixed, sc.lam_bar);
end


%  Step 2: the integrals and the cost aggregate. 
[Ia1, Ia2, Ix] = compute_integrals(w_min, w_min, sol_regime, sc);
[Sigma, Omega] = compute_Sigma_Omega(Ia1, Ia2, Ix, w_min, w_min, sc);

if Sigma <= 0
    warning('solve_joint_binding_bgp: Sigma <= 0 at w_min = %.6f. Returning failed.', w_min);
    sol = make_failed(w_min);
    return
end


%  Step 3: z = (psi_A*Omega^exp_Omega)^(1/exp_z)
z = (psi_A * Omega^sc.exp_Omega) ^ (1 / sc.exp_z);

if ~isreal(z) || z <= 0
    warning('solve_joint_binding_bgp: non-positive/complex z at w_min = %.6f. Returning failed.', w_min);
    sol = make_failed(w_min);
    return
end


%  Step 4: top layer
y   = sc.A * z^sc.exp_y;
k   = (sc.alpha / (sc.r + sc.delta)) * y;
w_n = (1 - sc.alpha - sc.psi) * y / sc.s_n;

if strcmp(sol_regime, 'no_auto')
    X = 0;
else
    X = (sc.psi * y / sc.p) * (sc.p^(-gam) * Ix / Sigma);
end

%  Step 5: labor markets
L_a1 = sc.psi * y * (w_min ^ (-1 - gam)) * Ia1 / Sigma;   % zero, the integral being zero
u_a1 = sc.s_a1 - L_a1;

L_a2 = sc.psi * y * (w_min ^ (-1 - gam)) * Ia2 / Sigma;
u_a2 = sc.s_a2 - L_a2;

% check
chk = struct();

chk.hat_eta_is_one = abs((w_min / w_min) - 1) < 1e-10;   
chk.eta_bar_gt_1    = sc.eta_bar > 1;                      
if strcmp(sol_regime, 'no_auto')
    chk.lambda_ok = true;
else
    chk.lambda_ok = (hat_lam_fixed < sc.lam_bar);
end
chk.wages_pos = (w_min > 0) && (w_n > 0);
chk.X_nonneg  = (X >= 0);
chk.L_a1_zero = abs(L_a1) < 1e-8;   

chk.core_ok = chk.hat_eta_is_one && chk.eta_bar_gt_1 && chk.lambda_ok && ...
              chk.wages_pos && chk.X_nonneg && chk.L_a1_zero;

if ~chk.core_ok
    warning('solve_joint_binding_bgp: core consistency check failed at w_min = %.6f. See sol.checks.', w_min);
    sol        = make_failed(w_min);
    sol.checks = chk;
    return
end

% Type-2 unemployment has to be non-negative
chk.u_a2_nonneg = (u_a2 >= 0);
if ~chk.u_a2_nonneg
    if params.verbose
        fprintf('FAILED: u_a2 = %.6f < 0. w_min not high enough to bind type-2. Returning failed.\n', u_a2);
    end
    sol        = make_failed(w_min);
    sol.checks = chk;
    return
end

chk.all_ok = chk.core_ok && chk.u_a2_nonneg;


%  output
sol.z          = z;
sol.w_a1       = w_min;
sol.w_a2       = w_min;
sol.hat_eta    = w_min / w_min;      % one, exactly
sol.hat_lambda = hat_lam_fixed;
sol.I_a1       = Ia1;
sol.I_a2       = Ia2;
sol.I_x        = Ix;
sol.Sigma      = Sigma;
sol.Omega      = Omega;
sol.y          = y;
sol.k          = k;
sol.w_n        = w_n;
sol.X          = X;
sol.L_a1       = L_a1;
sol.u_a1       = u_a1;
sol.L_a2       = L_a2;
sol.u_a2       = u_a2;
sol.regime     = sol_regime;
sol.w_min      = w_min;
sol.converged  = true;
sol.checks     = chk;

if params.verbose
    fprintf('\n--- Joint-Binding BGP Solution ---\n');
    fprintf('  Regime             : %s\n', sol_regime);
    fprintf('  w_min (=w_a1=w_a2) : %.6f\n', w_min);
    fprintf('  z                  : %.6f\n', z);
    fprintf('  y                  : %.6f\n', y);
    fprintf('  w_n                : %.6f\n', w_n);
    fprintf('  X (robots)         : %.6f\n', X);
    fprintf('  L_a2 = %.6f   u_a2 = %.6f\n', L_a2, u_a2);
    fprintf('  u_a1 = %.6f   (= s_a1, arithmetic)\n', u_a1);
    fprintf('\n');
end

end  % main function



function s = make_failed(w_min)
    s           = struct();
    s.converged = false;
    s.regime    = 'failed';
    s.w_min     = w_min;
    s.z = NaN; s.w_a1 = NaN; s.w_a2 = NaN;
    s.y = NaN; s.k = NaN; s.w_n = NaN; s.X = NaN;
    s.L_a1 = NaN; s.u_a1 = NaN; s.L_a2 = NaN; s.u_a2 = NaN;
    s.I_a1 = NaN; s.I_a2 = NaN; s.I_x = NaN; s.Sigma = NaN; s.Omega = NaN;
    s.hat_eta = NaN; s.hat_lambda = NaN;
end
