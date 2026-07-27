function value = residual_martensite_limit(p,state)
%SMA.RESIDUAL_MARTENSITE_LIMIT Cycle-evolved reverse-MT lower bound.

value = state.xi_res;
if ~p.enable_residual_martensite
    return;
end

peak_fraction = min(max(state.xi_r,0.0),1.0);
saturation = p.residual_xi_offset+ ...
    p.residual_xi_peak_gain*peak_fraction;
evolved = saturation*(1.0-exp(-p.residual_xi_a* ...
    max(state.zeta_d,0.0)));
value = min(max([state.xi_res,evolved,0.0]),state.xi);
end
