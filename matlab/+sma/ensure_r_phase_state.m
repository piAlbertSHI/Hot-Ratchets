function state = ensure_r_phase_state(p,state)
%SMA.ENSURE_R_PHASE_STATE Add and validate the committed R-phase fraction.

if ~isfield(state,'xi_R') || isempty(state.xi_R)
    state.xi_R=0.0;
end
assert(isscalar(state.xi_R) && isfinite(state.xi_R), ...
    'sma:ensure_r_phase_state:Finite', ...
    'The committed R-phase fraction must be a finite scalar.');
if p.enable_explicit_r_phase_state
    xi_M=min(max(state.xi,0.0),1.0);
    assert(state.xi_R>=-1.0e-12 && ...
        state.xi_R<=1.0-xi_M+1.0e-12, ...
        'sma:ensure_r_phase_state:Bounds', ...
        'The committed A/R/M phase fractions violate capacity.');
    state.xi_R=min(max(state.xi_R,0.0),1.0-xi_M);
end
end
