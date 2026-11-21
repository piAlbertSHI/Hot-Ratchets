function [xi_A,xi_R,xi_M,relative_R,branch,pi_R] = ...
    r_phase_return_map(p,xi_M,seq,T,training_activity, ...
    xi_M_old,xi_R_old)
%SMA.R_PHASE_RETURN_MAP Rate-independent Zhang Eq. (38) RT update.

xi_M=min(max(xi_M,0.0),1.0);
xi_M_old=min(max(xi_M_old,0.0),1.0);
old_inventory=max(1.0-xi_M_old,0.0);
if old_inventory>p.numerics.xi_tol
    relative_old=min(max(xi_R_old/old_inventory,0.0),1.0);
else
    relative_old=0.0;
end
training=1.0-exp(-p.R_training_rate*max(training_activity,0.0));
hardening=p.R_hardening/(1.0+p.R_training_gain*training);
base=p.R_H*max(seq,0.0)+p.R_entropy*(T-p.R_Tref)-p.R_tau0;
pi_R=base-hardening*relative_old;
tol=p.numerics.surface_tol;
branch=0.0;
relative_R=relative_old;
if pi_R>p.R_dissipation+tol && relative_old<1.0-p.numerics.xi_tol
    relative_R=(base-p.R_dissipation)/hardening;
    branch=1.0;
elseif pi_R<-p.R_dissipation-tol && ...
        relative_old>p.numerics.xi_tol
    relative_R=(base+p.R_dissipation)/hardening;
    branch=-1.0;
end
relative_R=min(max(relative_R,0.0),1.0);
xi_R=(1.0-xi_M)*relative_R;
xi_A=1.0-xi_M-xi_R;
xi_A=max(xi_A,0.0);
assert(abs(xi_A+xi_R+xi_M-1.0)<5.0e-14 && ...
    min([xi_A,xi_R,xi_M])>=-1.0e-14, ...
    'sma:r_phase_return_map:Conservation', ...
    'RT return map violates A/R/M phase conservation.');
end
