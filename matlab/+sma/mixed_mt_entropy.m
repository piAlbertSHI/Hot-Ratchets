function delta_eta_AR_M=mixed_mt_entropy(p,relative_R,delta_eta_p)
%SMA.MIXED_MT_ENTROPY Zhang Eq.(37-b), in MPa/K.
% eta_M - [(1-r_R)*eta_A+r_R*eta_R]
% = DeltaEta_A_M - r_R*DeltaEta_A_R.
% delta_eta_p is the independently owned configurational entropy shift.
if nargin<3, delta_eta_p=0; end
assert(all(isfinite(relative_R(:))) && all(isfinite(delta_eta_p(:))) && ...
    all(relative_R(:)>=-1e-12 & relative_R(:)<=1+1e-12), ...
    'sma:mixed_mt_entropy:Mixture','The A/R mixture weight must lie in [0,1].');
delta_eta_AR_M=p.rho_delta_s-relative_R.*p.R_entropy+delta_eta_p;
end
