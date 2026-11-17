function [g_AR_M,delta_eta_AR_M,relative_R,g_R_M] = ...
    mixed_mt_properties(p,xi_M,xi_R,g_A_M)
%SMA.MIXED_MT_PROPERTIES Zhang Eqs. (6) and (37-b).
% xi_R is the absolute R-phase fraction and g_A_M is the current A-to-M
% transformation-strain magnitude.  The identity g_A_R+g_R_M=g_A_M is
% enforced before applying the A/R mixture rule.

assert(isfinite(xi_M) && isfinite(xi_R) && isfinite(g_A_M), ...
    'Phase fractions and transformation strain must be finite.');
assert(xi_M>=-1.0e-12 && xi_R>=-1.0e-12 && ...
    xi_M+xi_R<=1.0+1.0e-12, ...
    'A/R/M phase fractions are not admissible.');
xi_M=min(max(xi_M,0.0),1.0);
inventory=max(1.0-xi_M,0.0);
xi_R=min(max(xi_R,0.0),inventory);
if inventory>p.numerics.xi_tol
    relative_R=xi_R/inventory;
else
    relative_R=0.0;
end
relative_R=min(max(relative_R,0.0),1.0);
g_R_M=max(g_A_M-p.R_H,0.0);
g_AR_M=(1.0-relative_R)*g_A_M+relative_R*g_R_M;
% At r_R=0 the parent is A; at r_R=1 it is R, so the entropy
% endpoints are DeltaEta_A_M and DeltaEta_R_M=DeltaEta_A_M-DeltaEta_A_R.
delta_eta_AR_M=sma.mixed_mt_entropy(p,relative_R);
end
