function value=candidate_zhang_source_mt_resistance( ...
    p,branch,xi,xi_recent,xi_res)
%SMA.CANDIDATE_ZHANG_SOURCE_MT_RESISTANCE Zhang Eqs. (61)-(63).
% This is a source-traceable shadow law; no observation or cycle index enters.
% History values are transaction inputs, not quantities to silently repair:
% clipping them can make an invalid residual-martensite owner look like a
% valid source state.  Permit only round-off at the physical bounds.
assert(isfield(p,'case_id') && p.case_id==1, ...
    'sma:candidate_zhang_source_mt_resistance:MaterialRole', ...
    'Zhang source MT resistance is defined only for the NiTi Zhang-2022 card.');
tol=1e-12;
assert(isfinite(branch) && isfinite(xi) && isfinite(xi_recent) && ...
    isfinite(xi_res),'Zhang source MT inputs must be finite.');
assert(branch==1 || branch==-1,'Zhang source MT branch must be +/-1.');
assert(xi>=-tol && xi<=1+tol,'Current M fraction is outside [0,1].');
assert(xi_res>=-tol && xi_res<=xi+tol, ...
    'Residual M fraction must satisfy 0<=xi_res<=xi.');
xi=min(max(xi,0.0),1.0);
xi_res=min(max(xi_res,0.0),xi);
assert(xi_recent>=-tol && xi_recent<=1+tol, ...
    'xi_M_recent is outside the physical source-history bounds.');
xi_recent=min(max(xi_recent,0.0),1.0);
if branch<0
    % Eq. (62) reverse hardening divides by the total M fraction at the
    % beginning of the reverse episode. That episode memory must dominate
    % both the residual-M owner and the current total M. The forward branch
    % does not use xi_recent in tau_s; its xi_eff remains well-defined even
    % when Eq. (53-b) generates residual M in the same forward episode.
    assert(xi_recent>=xi_res-tol, ...
        'Reverse MT requires xi_M_recent >= xi_Mr.');
    assert(xi_recent>=xi-tol, ...
        'Reverse MT requires xi_M_recent >= current xi_M.');
end
denominator=max(1.0-xi_res,eps);
xi_eff=min(max((xi-xi_res)/denominator,0.0),1.0);
if branch>0
    nonlinear=p.zhang_hfor_AR_M*xi_eff^p.zhang_ntr_AR_M;
else
    % Zhang Eq. (62): reverse hardening depends on the current effective
    % M fraction relative to the M fraction at the beginning of reverse
    % transformation.  It is not a constant based only on xi_recent.
    reverse_ratio=min(max(xi_eff/max(xi_recent,eps),0.0),1.0);
    nonlinear=p.zhang_hrev_AR_M*reverse_ratio^p.zhang_ntr_AR_M;
end
value=p.zhang_tau0_AR_M+p.zhang_H_AR_M*xi_eff+nonlinear;
end
