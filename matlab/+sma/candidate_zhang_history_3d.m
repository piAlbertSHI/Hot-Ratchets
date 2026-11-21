function history = candidate_zhang_history_3d(p,state,xi,h_tp,tau)
%SMA.CANDIDATE_ZHANG_HISTORY_3D Candidate-coupled tensor history.

history=sma.zhang_direct_history_from_state(state);
if p.zhang_direct_Yp_coupling<=0 && ...
        p.zhang_direct_Omega_coupling<=0, return; end
h_arp=sma.update_arp_3d(p,state,tau,xi);
darp=h_arp-state.h_arp; dtp=h_tp-state.h_tp;
gamma_arp=sqrt(max((2.0/3.0)*sum(darp.^2,'all'),0.0));
gamma_tp=sqrt(max((2.0/3.0)*sum(dtp.^2,'all'),0.0));
history=sma.zhang_history_evolution_direct( ...
    history,xi,gamma_arp,gamma_tp);
end
