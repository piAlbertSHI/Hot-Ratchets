function [phase,increment]=candidate_zhang_current_mixture_phase_increment_1d(p,old,dR,dM)
%SMA.CANDIDATE_ZHANG_CURRENT_MIXTURE_PHASE_INCREMENT_1D
% Endpoint-weight backward Euler discretization of Zhang Eqs. (4)--(7)
% with the current A/(A+R) mixture. It is consistent with the rate law,
% but is not the exact finite-step solution for constant signed rates.
% xi_MA/xi_MR remain passive origin labels; xi_M is primary.

assert(all(isfinite([old.xi_A,old.xi_R,old.xi_M,dR,dM])), ...
    'sma:current_mixture:Finite','Phase state and increments must be finite.');
mixed=old.xi_A+old.xi_R;
assert(mixed>0,'sma:current_mixture:UndefinedPureM', ...
    'Current A/(A+R) composition is undefined at pure M.');
[phase,increment]=sma.candidate_coupled_phase_increment(p,old,dR,dM);
wR=(old.xi_R+dR)/mixed;
wA=(old.xi_A-dR)/mixed;
assert(wA>=-50*eps && wR>=-50*eps && abs(wA+wR-1)<=5e-14, ...
    'sma:current_mixture:Composition','Current A/R mixture is inadmissible.');
if dM<0
    assert(old.xi_M>0 && dM>=-old.xi_M, ...
        'sma:current_mixture:ReverseCapacity','Reverse flow exceeds M inventory.');
    fraction=1+dM/old.xi_M;
    phase.xi_MA=old.xi_MA*fraction;
    phase.xi_MR=old.xi_MR*fraction;
    phase.xi_M=old.xi_M+dM;
    phase.xi_A=old.xi_A-dR-wA*dM;
    phase.xi_R=old.xi_R+dR-wR*dM;
    increment.delta_xi_MA=phase.xi_MA-old.xi_MA;
    increment.delta_xi_MR=phase.xi_MR-old.xi_MR;
    increment.weight_A=wA;
    increment.weight_R=wR;
    increment.g_mixed=wA*p.Hmax_i+wR*(p.Hmax_i-p.R_H);
    increment.delta_eps_MT=increment.g_mixed*dM;
    phase.eps_MT=old.eps_MT+increment.delta_eps_MT;
end
phase.xi_M=old.xi_M+dM;
increment.delta_lambda_AM=wA*dM;
increment.delta_lambda_RM=wR*dM;
increment.composition_mode='source_current_mixture';
increment.passive_origin_labels=true;
increment.passive_origin_sum_defect=phase.xi_MA+phase.xi_MR-phase.xi_M;
increment.phase_conservation_error=phase.xi_A+phase.xi_R+phase.xi_M-1;
increment.source_be_residual=[phase.xi_A-old.xi_A+dR+wA*dM; ...
    phase.xi_R-old.xi_R-dR+wR*dM;phase.xi_M-old.xi_M-dM];
increment.strain_inventory_defect=increment.delta_eps_RT+increment.delta_eps_MT- ...
    p.R_H*(phase.xi_R-old.xi_R)-p.Hmax_i*(phase.xi_M-old.xi_M);
assert(max(abs(increment.source_be_residual))<5e-14 && ...
    abs(increment.phase_conservation_error)<5e-14 && ...
    abs(increment.strain_inventory_defect)<5e-14 && ...
    min([phase.xi_A,phase.xi_R,phase.xi_M,phase.xi_MA,phase.xi_MR])>=-50*eps, ...
    'sma:current_mixture:Conservation','Current-mixture transfer failed.');
end
