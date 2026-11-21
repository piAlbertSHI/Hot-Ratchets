function [phi_f,phi_r] = trial_surfaces_1d(p,state,sigma,T)
%SMA.TRIAL_SURFACES_1D Evaluate both transformation branches at trial state.

sigma_eff = sigma+state.beta;
if abs(sigma_eff)>p.numerics.direction_tol
    direction = sign(sigma_eff);
else
    direction = state.direction;
end
[Hcur,Hmax] = sma.hcur(p,abs(sigma_eff),state.zeta_d);
ratio = Hcur/max(Hmax,eps);
[Hcur_mixed,delta_eta_AR_M]=sma.candidate_mixed_mt_properties( ...
    p,state,state.xi,abs(sigma),T,Hcur);
gp = p.Cp1*p.Cp2/(1.0+p.Cp2*state.zeta_d);
lambda_f = Hcur_mixed*direction;
lambda_tp_f = ratio^2*gp*direction;
[~,phi_f] = sma.driving_force_1d( ...
    p,1,sigma,state.beta,T,state.xi,lambda_f,lambda_tp_f, ...
    state.zeta_d,state.xi_r,state.zeta_p,state.xi_res, ...
    state.Y_population_activity,state.zhang_Yp,state.zhang_Omega, ...
    delta_eta_AR_M);
if state.xi_r>p.numerics.xi_tol
    lambda_r = state.eps_tr_r/state.xi_r;
else
    lambda_r = 0.0;
end
lambda_tp_r = -ratio^2*gp*lambda_r;
[~,phi_r] = sma.driving_force_1d( ...
    p,-1,sigma,state.beta,T,state.xi,lambda_r,lambda_tp_r, ...
    state.zeta_d,state.xi_r,state.zeta_p,state.xi_res, ...
    state.Y_population_activity,state.zhang_Yp,state.zhang_Omega, ...
    delta_eta_AR_M);
end
