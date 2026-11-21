function x = state_vector_1d_v2(s)
%SMA.STATE_VECTOR_1D_V2 Pack base state and threshold-population history.

if ~isfield(s,'Y_population_activity')
    error('sma:state_vector_1d_v2:MissingPopulationState', ...
        'The v2 state requires Y_population_activity.');
end
activity=s.Y_population_activity(:);
if any(~isfinite(activity)) || any(activity<0)
    error('sma:state_vector_1d_v2:InvalidPopulationState', ...
        'Population activities must be finite and nonnegative.');
end
% Construct the legacy prefix directly because state_vector_1d correctly
% rejects a lossy conversion when the population history is nonzero.
base=[s.T; s.xi; s.eps_tr; s.eps_tp; s.zeta_d; s.beta; ...
    s.sigma; s.eps; s.branch; s.direction; s.xi_r; s.eps_tr_r; ...
    s.eps_arp; s.xi_res; s.zeta_p];
x=[base;activity];
end
