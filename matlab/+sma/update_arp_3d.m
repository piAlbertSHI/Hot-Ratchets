function value = update_arp_3d(p,state,tau,xi)
%SMA.UPDATE_ARP_3D J2 A/R-phase plastic correction.

value = state.h_arp;
if ~p.enable_arp
    return;
end

[seq,dev] = sma.equivalent_stress(tau);
if seq<=p.numerics.direction_tol
    return;
end
gamma_old = sqrt(2.0/3.0)*norm(state.h_arp,'fro');
peak_fraction = min(max(max(state.xi_r,xi),0.0),1.0);
tau_progress = min(max((peak_fraction-p.arp_tau_fraction_start)/ ...
    p.arp_tau_fraction_width,0.0),1.0);
tau_progress = tau_progress^2*(3.0-2.0*tau_progress);
tau0 = p.arp_tau0+(p.arp_tau0_f-p.arp_tau0)*tau_progress;
gamma_surface = max((seq-tau0)/p.arp_H,0.0);
matrix_fraction = max(1.0-min(max(xi,0.0),1.0),0.0)^ ...
    p.arp_matrix_power;
delta_gamma = matrix_fraction*max(gamma_surface-gamma_old,0.0);
N = 1.5*dev/seq;
value = state.h_arp+delta_gamma*N;
end
