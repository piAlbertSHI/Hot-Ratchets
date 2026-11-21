function state = unpack_state_1d(p,x)
%SMA.UNPACK_STATE_1D Unpack all history required by active mechanisms.

x=x(:);
assert(numel(x)==sma.state_width_1d(p), ...
    'sma:unpack_state_1d:Width','State width does not match parameters.');
state=sma.state_from_vector_1d(x(1:15),p);
cursor=16;
if p.Y_activity_population_mix>0
    n=p.Y_activity_population_bins;
    state.Y_population_activity=x(cursor:cursor+n-1); cursor=cursor+n;
end
if p.enable_explicit_r_phase_state
    state.xi_R=x(cursor); cursor=cursor+1;
    state=sma.ensure_r_phase_state(p,state);
end
if p.zhang_direct_Yp_coupling>0 || p.zhang_direct_Omega_coupling>0
    z=x(cursor:cursor+8);
    names={'zhang_xi_Mr','zhang_rho_M','zhang_rho_AR','zhang_D_M', ...
        'zhang_D_AR','zhang_DeltaEtaP','zhang_Omega','zhang_Yp', ...
        'zhang_Lambda_M'};
    for k=1:numel(names), state.(names{k})=z(k); end
    cursor=cursor+9;
end
assert(cursor==numel(x)+1,'sma:unpack_state_1d:Cursor', ...
    'State unpack cursor mismatch.');
end
