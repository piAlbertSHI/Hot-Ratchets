function [phi_f,phi_r] = trial_surfaces_3d(p,state,tau,T)
%SMA.TRIAL_SURFACES_3D Evaluate both tensor transformation branches.

[seq,dev] = sma.equivalent_stress(tau+state.beta);
if seq>p.numerics.direction_tol
    N = 1.5*dev/seq;
else
    N = state.N;
end
[Hcur,Hmax] = sma.hcur(p,seq,state.zeta_d);
ratio = Hcur/max(Hmax,eps);
[seq_R,~]=sma.equivalent_stress(tau);
[Hcur_mixed,delta_eta_AR_M]=sma.candidate_mixed_mt_properties( ...
    p,state,state.xi,seq_R,T,Hcur);
gp = p.Cp1*p.Cp2/(1.0+p.Cp2*state.zeta_d);
Lambda_f = Hcur_mixed*N;
Lambda_tp_f = ratio^2*gp*N;
[~,phi_f] = sma.driving_force_3d( ...
    p,1,tau,state.beta,T,state.xi,Lambda_f,Lambda_tp_f, ...
    state.zeta_d,state.xi_r,state.zeta_p,state.xi_res, ...
    state.Y_population_activity,state.zhang_Yp,state.zhang_Omega, ...
    delta_eta_AR_M);
if state.xi_r>p.numerics.xi_tol
    Lambda_r = state.h_tr_r/state.xi_r;
else
    Lambda_r = zeros(3);
end
Lambda_tp_r = -ratio^2*gp*Lambda_r;
[~,phi_r] = sma.driving_force_3d( ...
    p,-1,tau,state.beta,T,state.xi,Lambda_r,Lambda_tp_r, ...
    state.zeta_d,state.xi_r,state.zeta_p,state.xi_res, ...
    state.Y_population_activity,state.zhang_Yp,state.zhang_Omega, ...
    delta_eta_AR_M);
end
