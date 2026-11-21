function state = ensure_zhang_direct_state(p,state)
%SMA.ENSURE_ZHANG_DIRECT_STATE Add/validate versioned direct history.

names={'zhang_xi_Mr','zhang_rho_M','zhang_rho_AR','zhang_D_M', ...
    'zhang_D_AR','zhang_DeltaEtaP','zhang_Omega','zhang_Yp', ...
    'zhang_Lambda_M'};
missing=~all(isfield(state,names));
if missing && (p.zhang_direct_Yp_coupling>0 || ...
        p.zhang_direct_Omega_coupling>0)
    error('sma:ensure_zhang_direct_state:MissingActiveState', ...
        'Enabled direct Zhang coupling requires a versioned state.');
end
for k=1:numel(names)
    if ~isfield(state,names{k}), state.(names{k})=0.0; end
end
v=cellfun(@(n)state.(n),names);
if any(~isfinite(v)) || any(v([1:5,9])<0)
    error('sma:ensure_zhang_direct_state:InvalidState', ...
        'Direct Zhang history is invalid.');
end
end
