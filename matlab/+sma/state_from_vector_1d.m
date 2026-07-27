function s = state_from_vector_1d(x, p)
%SMA.STATE_FROM_VECTOR_1D Unpack the standard 15-value Simulink state.

x = x(:);
if numel(x) ~= 15
    error('sma:state_from_vector_1d:Width', ...
        'The standard 1D state must contain 15 entries.');
end
s = sma.initial_state_1d(p);
names = {'T','xi','eps_tr','eps_tp','zeta_d','beta','sigma','eps', ...
    'branch','direction','xi_r','eps_tr_r','eps_arp','xi_res','zeta_p'};
for k = 1:numel(names)
    s.(names{k}) = x(k);
end
[s.Hcur,s.Hmax] = sma.hcur(p,abs(s.sigma+s.beta),s.zeta_d);
s.xi_d = min(max(s.Hcur/max(s.Hmax,eps),0.0),1.0)*s.xi;
end
