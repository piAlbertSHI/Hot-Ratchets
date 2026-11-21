function [R, data] = residual_1d(p, state, sigma, T, branch, reverse, x)
%SMA.RESIDUAL_1D Coupled local residual for the strict 1D reduction.
% x = [xi, eps_tr, eps_tp, zeta_d, beta].

xi = x(1);
eps_tr = x(2);
eps_tp = x(3);
zeta_d = x(4);
beta = x(5);
sigma_eff = sigma + beta;
if branch < 0 && abs(reverse.eps_tr_r) > p.numerics.direction_tol
    direction = sign(reverse.eps_tr_r);
elseif abs(sigma_eff) > p.numerics.direction_tol
    direction = sign(sigma_eff);
else
    direction = state.direction;
end
if direction == 0.0
    direction = 1.0;
end

[Hcur, Hmax] = sma.hcur(p,abs(sigma_eff),zeta_d);
ratio = min(max(Hcur / max(Hmax, eps), 0.0), 1.0);
[Hcur_mixed,delta_eta_AR_M,relative_R_for_MT,xi_R_for_MT] = ...
    sma.candidate_mixed_mt_properties( ...
    p,state,xi,abs(sigma),T,Hcur);
if branch > 0
    lambda_tr = Hcur_mixed * direction;
else
    if reverse.xi_r > p.numerics.xi_tol
        lambda_tr = reverse.eps_tr_r / reverse.xi_r;
    else
        lambda_tr = 0.0;
    end
end

gp = p.Cp1 * p.Cp2 / (1.0 + p.Cp2 * max(zeta_d, 0.0));
if branch > 0
    lambda_tp = ratio^2 * gp * direction;
else
    if reverse.xi_r > p.numerics.xi_tol
        lambda_tp = -ratio^2 * gp * reverse.eps_tr_r / reverse.xi_r;
    else
        lambda_tp = 0.0;
    end
end

xi_d = ratio * xi;
beta_target = p.sigma_b * ratio * ...
    (1.0 - exp(-p.lambda_beta * max(zeta_d, 0.0))) * direction;
population_activity=sma.update_population_activity( ...
    p,state.Y_population_activity,state.xi,xi);
zhang_history=sma.candidate_zhang_history_1d( ...
    p,state,xi,eps_tp,sigma);
[pi_value, phi, force_terms] = sma.driving_force_1d( ...
    p,branch,sigma,beta,T,xi,lambda_tr,lambda_tp,zeta_d, ...
    reverse.xi_r,state.zeta_p,state.xi_res,population_activity, ...
    zhang_history.Yp,zhang_history.Omega,delta_eta_AR_M);
threshold = sma.transformation_resistance( ...
    p,branch,zeta_d,xi,reverse.xi_r,state.zeta_p,state.xi_res, ...
    population_activity,zhang_history.Yp);
delta_xi = xi - state.xi;
delta_xi_d = xi_d - state.xi_d;

strain_scale = max(p.Hmax_i, 1.0e-4);
trip_scale = max(p.Cp1, 1.0e-4);
beta_scale = max([p.sigma_b, p.reference_stress, 1.0]);
R = [
    phi / max(p.reference_stress, 1.0)
    (eps_tr - state.eps_tr - lambda_tr * delta_xi) / strain_scale
    (eps_tp - state.eps_tp - lambda_tp * delta_xi) / trip_scale
    zeta_d - state.zeta_d - abs(delta_xi_d)
    (beta - beta_target) / beta_scale
    ];

data = struct( ...
    'phi', phi, ...
    'pi', pi_value, ...
    'force_terms', force_terms, ...
    'threshold', threshold, ...
    'direction', direction, ...
    'Hcur', Hcur, ...
    'Hcur_mixed', Hcur_mixed, ...
    'Hmax', Hmax, ...
    'delta_eta_AR_M',delta_eta_AR_M, ...
    'relative_R_for_MT',relative_R_for_MT, ...
    'xi_R_for_MT',xi_R_for_MT, ...
    'ratio', ratio, ...
    'lambda_tr', lambda_tr, ...
    'lambda_tp', lambda_tp, ...
    'xi_d', xi_d, ...
    'Y_population_activity',population_activity, ...
    'zhang_history',zhang_history, ...
    'beta_target', beta_target);
end
