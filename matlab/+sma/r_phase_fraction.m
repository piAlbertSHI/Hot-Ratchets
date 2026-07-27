function [xi_R,dxi_dseq] = r_phase_fraction( ...
    p,seq,T,training_activity)
%SMA.R_PHASE_FRACTION Reduced reversible A-to-R equilibrium branch.

if nargin<4
    training_activity = 0.0;
end
if ~p.enable_r_phase || p.R_H<=0.0
    xi_R = 0.0;
    dxi_dseq = 0.0;
    return;
end

training = 1.0-exp(-p.R_training_rate*max(training_activity,0.0));
hardening = p.R_hardening/(1.0+p.R_training_gain*training);
driving = p.R_H*max(seq,0.0)+ ...
    p.R_entropy*(T-p.R_Tref)-p.R_tau0;
raw_fraction = driving/hardening;
xi_R = min(max(raw_fraction,0.0),p.R_max_fraction);
if raw_fraction>0.0 && raw_fraction<p.R_max_fraction
    dxi_dseq = p.R_H/hardening;
else
    dxi_dseq = 0.0;
end
end
