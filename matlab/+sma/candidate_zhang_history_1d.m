function history = candidate_zhang_history_1d(p,state,xi,eps_tp,sigma)
%SMA.CANDIDATE_ZHANG_HISTORY_1D Candidate-coupled direct history.

history=sma.zhang_direct_history_from_state(state);
if p.zhang_direct_Yp_coupling<=0 && ...
        p.zhang_direct_Omega_coupling<=0, return; end
eps_arp=sma.update_arp_1d(p,state,sigma,xi);
history=sma.zhang_history_evolution_direct(history,xi, ...
    eps_arp-state.eps_arp,eps_tp-state.eps_tp);
end
