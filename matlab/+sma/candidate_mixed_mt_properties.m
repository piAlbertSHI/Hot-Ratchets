function [g_AR_M,delta_eta_AR_M,relative_R,xi_R] = ...
    candidate_mixed_mt_properties(p,state,xi_M,seq,T,g_A_M)
%SMA.CANDIDATE_MIXED_MT_PROPERTIES Trial-state Zhang Eqs. (6)/(37-b).
% The independent RT map is analytically eliminated at each MT residual
% probe.  No observation, cycle number or material-case branch is read.

if ~p.enable_composition_weighted_mt
    g_AR_M=g_A_M;
    delta_eta_AR_M=p.rho_delta_s;
    relative_R=0.0;
    xi_R=state.xi_R;
    return;
end
assert(p.enable_explicit_r_phase_state && p.enable_r_phase_return_map, ...
    'sma:candidate_mixed_mt_properties:ExplicitRTRequired', ...
    'Composition-weighted MT requires the independent explicit RT map.');
state=sma.ensure_r_phase_state(p,state);
[~,xi_R,xi_M,relative_map]=sma.r_phase_return_map( ...
    p,xi_M,max(seq,0.0),T,state.zeta_p,state.xi,state.xi_R);
[g_AR_M,delta_eta_AR_M,relative_check]= ...
    sma.mixed_mt_properties(p,xi_M,xi_R,g_A_M);
if 1.0-xi_M>p.numerics.xi_tol
    assert(abs(relative_map-relative_check)<5.0e-12, ...
        'RT elimination and mixed-MT inventory are inconsistent.');
end
relative_R=relative_check;
end
