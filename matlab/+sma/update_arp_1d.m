function value = update_arp_1d(p,state,sigma,xi)
%SMA.UPDATE_ARP_1D Rate-independent A/R-phase plastic correction.

value = state.eps_arp;
if ~p.enable_arp
    return;
end

direction = sign(sigma);
if direction==0.0
    direction = sign(state.eps_arp);
end
if direction==0.0
    direction = 1.0;
end

gamma_old = abs(state.eps_arp);
peak_fraction = min(max(max(state.xi_r,xi),0.0),1.0);
tau_progress = min(max((peak_fraction-p.arp_tau_fraction_start)/ ...
    p.arp_tau_fraction_width,0.0),1.0);
tau_progress = tau_progress^2*(3.0-2.0*tau_progress);
tau0 = p.arp_tau0+(p.arp_tau0_f-p.arp_tau0)*tau_progress;
gamma_surface = max((abs(sigma)-tau0)/p.arp_H,0.0);
matrix_fraction = max(1.0-min(max(xi,0.0),1.0),0.0)^ ...
    p.arp_matrix_power;
delta_gamma = matrix_fraction*max(gamma_surface-gamma_old,0.0);
value = state.eps_arp+direction*delta_gamma;
end
