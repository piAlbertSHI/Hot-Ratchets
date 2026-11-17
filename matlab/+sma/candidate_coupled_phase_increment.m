function [phase,increment] = candidate_coupled_phase_increment( ...
    p,old,delta_lambda_R,delta_lambda_M)
%SMA.CANDIDATE_COUPLED_PHASE_INCREMENT Origin-memory inventory step.
% This pure shadow kernel has no observation access and is not called by the
% production return map. RT acts first. Forward MT uses that A/R mixture;
% reverse MT returns the stored martensite-origin proportions. The latter
% is a project extension and need not satisfy Zhang Eq. (7)'s instantaneous
% current A/R composition when the stored and current ratios differ.

required={'xi_A','xi_R','xi_M','xi_MA','xi_MR','eps_RT','eps_MT'};
for k=1:numel(required)
    assert(isfield(old,required{k}), ...
        'sma:candidate_coupled_phase_increment:State', ...
        'Missing phase state field %s.',required{k});
end
values=[old.xi_A,old.xi_R,old.xi_M,old.xi_MA,old.xi_MR, ...
    old.eps_RT,old.eps_MT, ...
    delta_lambda_R,delta_lambda_M];
assert(all(isfinite(values)), ...
    'sma:candidate_coupled_phase_increment:Finite', ...
    'Phase state and increments must be finite.');
tol=50.0*eps;
phase_values=[old.xi_A,old.xi_R,old.xi_M,old.xi_MA,old.xi_MR];
old_phase_error_raw=old.xi_A+old.xi_R+old.xi_M-1.0;
old_source_error_raw=old.xi_MA+old.xi_MR-old.xi_M;
% A committed state may cross the strict 50*eps gate after many accepted
% transfers even though the underlying inventory is still physical.  Repair
% only this representational roundoff: the M source labels own xi_M and A is
% reconstructed as the exact inventory complement.  Any larger discrepancy
% remains a hard admissibility failure.
roundoff_tol=256.0*eps*max([1.0,abs(phase_values)]);
old_roundoff_correction=[0.0,0.0,0.0];
old_roundoff_canonicalized=false;
if min(phase_values)>=-tol && ...
        abs(old_phase_error_raw)<=roundoff_tol && ...
        abs(old_source_error_raw)<=roundoff_tol && ...
        (old_phase_error_raw~=0.0 || old_source_error_raw~=0.0)
    old_xi_A=old.xi_A;
    old_xi_M=old.xi_M;
    old.xi_M=old.xi_MA+old.xi_MR;
    old.xi_A=1.0-old.xi_R-old.xi_M;
    old_roundoff_correction=[old.xi_A-old_xi_A,old.xi_M-old_xi_M, ...
        old_source_error_raw];
    assert(max(abs(old_roundoff_correction))<=roundoff_tol, ...
        'sma:candidate_coupled_phase_increment:RoundoffCanonicalization', ...
        'Roundoff canonicalization exceeded its explicit bound.');
    old_roundoff_canonicalized=true;
else
    assert(min(phase_values)>=-tol && ...
        abs(old_phase_error_raw)<=tol && ...
        abs(old_source_error_raw)<=tol, ...
        'sma:candidate_coupled_phase_increment:OldConservation', ...
        'Committed A/R/M inventory is not admissible.');
end

% RT transfer at fixed mixed A+R inventory.
xi_A_after_RT=old.xi_A-delta_lambda_R;
xi_R_after_RT=old.xi_R+delta_lambda_R;
mixed_old=old.xi_A+old.xi_R;
assert(min([xi_A_after_RT,xi_R_after_RT])>=-tol, ...
    'sma:candidate_coupled_phase_increment:RTCapacity', ...
    'RT increment exceeds available A or R inventory.');

if abs(delta_lambda_M)<=tol
    if mixed_old>tol
        weight_A=xi_A_after_RT/mixed_old;
        weight_R=xi_R_after_RT/mixed_old;
    elseif old.xi_M>tol
        weight_A=old.xi_MA/old.xi_M;
        weight_R=old.xi_MR/old.xi_M;
    else
        weight_A=1.0;
        weight_R=0.0;
    end
elseif delta_lambda_M>0.0
    assert(delta_lambda_M<=mixed_old+tol && mixed_old>tol, ...
        'sma:candidate_coupled_phase_increment:ForwardMTCapacity', ...
        'Forward MT increment exceeds mixed A+R inventory.');
    weight_A=xi_A_after_RT/mixed_old;
    weight_R=xi_R_after_RT/mixed_old;
else
    assert(-delta_lambda_M<=old.xi_M+tol && old.xi_M>tol, ...
        'sma:candidate_coupled_phase_increment:ReverseMTCapacity', ...
        'Reverse MT increment exceeds M inventory.');
    % Reverse MT must return the stored source composition of martensite;
    % the remaining A/R mixture alone cannot identify this information.
    weight_A=old.xi_MA/old.xi_M;
    weight_R=old.xi_MR/old.xi_M;
end
delta_xi_MA=weight_A*delta_lambda_M;
delta_xi_MR=weight_R*delta_lambda_M;
xi_MA=old.xi_MA+delta_xi_MA;
xi_MR=old.xi_MR+delta_xi_MR;
xi_A=xi_A_after_RT-delta_xi_MA;
xi_R=xi_R_after_RT-delta_xi_MR;

xi_A=max(xi_A,0.0); xi_R=max(xi_R,0.0);
xi_MA=max(xi_MA,0.0); xi_MR=max(xi_MR,0.0);
xi_M_raw=xi_MA+xi_MR;
raw_phase_error=xi_A+xi_R+xi_M_raw-1.0;
final_roundoff_tol=256.0*eps*max([1.0,abs([xi_A,xi_R,xi_MA,xi_MR])]);
assert(abs(raw_phase_error)<=final_roundoff_tol, ...
    'sma:candidate_coupled_phase_increment:Conservation', ...
    'Accepted A/R/M increment requires more than roundoff canonicalization.');
xi_M=xi_M_raw;
xi_A_raw=xi_A;
xi_A=1.0-xi_R-xi_M;
final_roundoff_correction=xi_A-xi_A_raw;
assert(abs(final_roundoff_correction)<=final_roundoff_tol, ...
    'sma:candidate_coupled_phase_increment:Conservation', ...
    'Final inventory canonicalization exceeded its explicit bound.');
g_A_M=p.Hmax_i;
g_A_R=p.R_H;
g_R_M=g_A_M-g_A_R;
assert(g_A_M>=g_A_R && g_A_R>=0.0, ...
    'sma:candidate_coupled_phase_increment:TransformationStrain', ...
    'The source transformation-strain endpoints are not admissible.');
g_mixed=weight_A*g_A_M+weight_R*g_R_M;
delta_eps_RT=g_A_R*delta_lambda_R;
delta_eps_MT=g_A_M*delta_xi_MA+g_R_M*delta_xi_MR;

phase=old;
phase.xi_A=xi_A;
phase.xi_R=xi_R;
phase.xi_M=xi_M;
phase.xi_MA=xi_MA;
phase.xi_MR=xi_MR;
phase.eps_RT=old.eps_RT+delta_eps_RT;
phase.eps_MT=old.eps_MT+delta_eps_MT;
increment=struct('delta_lambda_R',delta_lambda_R, ...
    'delta_lambda_M',delta_lambda_M,'weight_A',weight_A, ...
    'weight_R',weight_R,'g_A_R',g_A_R,'g_R_M',g_R_M, ...
    'g_mixed',g_mixed,'delta_eps_RT',delta_eps_RT, ...
    'delta_eps_MT',delta_eps_MT,'delta_xi_MA',delta_xi_MA, ...
    'delta_xi_MR',delta_xi_MR, ...
    'old_phase_conservation_error_raw',old_phase_error_raw, ...
    'old_source_conservation_error_raw',old_source_error_raw, ...
    'old_roundoff_canonicalized',old_roundoff_canonicalized, ...
    'old_roundoff_correction',old_roundoff_correction, ...
    'final_roundoff_correction',final_roundoff_correction, ...
    'phase_conservation_error',xi_A+xi_R+xi_M-1.0);
assert(abs(increment.phase_conservation_error)<=5.0e-14 && ...
    min([xi_A,xi_R,xi_M,xi_MA,xi_MR])>=-tol && ...
    abs(xi_MA+xi_MR-xi_M)<=5.0e-14, ...
    'sma:candidate_coupled_phase_increment:Conservation', ...
    'Accepted A/R/M increment violates phase conservation.');
end
