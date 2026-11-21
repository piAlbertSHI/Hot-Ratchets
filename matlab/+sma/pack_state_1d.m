function x = pack_state_1d(p,state)
%SMA.PACK_STATE_1D Pack all history required by active mechanisms.

base=[state.T;state.xi;state.eps_tr;state.eps_tp;state.zeta_d; ...
    state.beta;state.sigma;state.eps;state.branch;state.direction; ...
    state.xi_r;state.eps_tr_r;state.eps_arp;state.xi_res;state.zeta_p];
x=base;
if p.Y_activity_population_mix>0
    activity=state.Y_population_activity(:);
    assert(numel(activity)==p.Y_activity_population_bins && ...
        all(isfinite(activity)) && all(activity>=0), ...
        'sma:pack_state_1d:Population','Invalid population history.');
    x=[x;activity]; %#ok<AGROW>
end
if p.enable_explicit_r_phase_state
    state=sma.ensure_r_phase_state(p,state);
    x=[x;state.xi_R]; %#ok<AGROW>
end
if p.zhang_direct_Yp_coupling>0 || p.zhang_direct_Omega_coupling>0
    z=[state.zhang_xi_Mr;state.zhang_rho_M;state.zhang_rho_AR; ...
        state.zhang_D_M;state.zhang_D_AR;state.zhang_DeltaEtaP; ...
        state.zhang_Omega;state.zhang_Yp;state.zhang_Lambda_M];
    assert(all(isfinite(z)) && all(z([1:5,9])>=0), ...
        'sma:pack_state_1d:Zhang','Invalid direct Zhang history.');
    x=[x;z]; %#ok<AGROW>
end
assert(numel(x)==sma.state_width_1d(p), ...
    'sma:pack_state_1d:Width','Packed state width mismatch.');
end
