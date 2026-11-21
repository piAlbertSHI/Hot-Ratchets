function width = state_width_1d(p)
%SMA.STATE_WIDTH_1D Versioned state width selected by active mechanisms.

width=15;
if p.Y_activity_population_mix>0
    width=width+p.Y_activity_population_bins;
end
if p.enable_explicit_r_phase_state
    width=width+1;
end
if p.zhang_direct_Yp_coupling>0 || p.zhang_direct_Omega_coupling>0
    width=width+9;
end
end
