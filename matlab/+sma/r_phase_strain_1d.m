function [eps_R,xi_R,deps_dsigma] = r_phase_strain_1d( ...
    p,sigma,T,training_activity,xi_M,xi_R_old,xi_M_old)
%SMA.R_PHASE_STRAIN_1D Reversible uniaxial R-phase strain.

if nargin<4
    training_activity = 0.0;
end
if nargin<5
    xi_M=0.0;
end
[relative_R,dxi_dseq] = sma.r_phase_fraction( ...
    p,abs(sigma),T,training_activity);
if p.enable_r_phase_return_map
    assert(nargin>=7,'sma:r_phase_strain_1d:CommittedState', ...
        'The RT return map requires committed R and M fractions.');
    [~,xi_R,xi_M,~,branch] = sma.r_phase_return_map( ...
        p,xi_M,abs(sigma),T,training_activity,xi_M_old,xi_R_old); %#ok<ASGLU>
    if branch==0
        dxi_dseq=0.0;
    else
        training=1.0-exp(-p.R_training_rate*max(training_activity,0.0));
        hardening=p.R_hardening/(1.0+p.R_training_gain*training);
        dxi_dseq=(1.0-xi_M)*p.R_H/hardening;
    end
elseif p.enable_explicit_r_phase_state
    [~,xi_R,xi_M] = sma.explicit_phase_fractions( ...
        p,xi_M,abs(sigma),T,training_activity); %#ok<ASGLU>
    dxi_dseq=(1.0-xi_M)*dxi_dseq;
else
    % Even on the legacy algebraic branch, the mechanical strain is a
    % mixture quantity and must use the absolute R inventory.  The previous
    % code returned relative_R here while r_phase_latent_specific multiplied
    % by (1-xi_M), so mechanics and latent heat used different phase
    % measures.  Keep the output xi_R in the same absolute convention as the
    % explicit branch and apply the corresponding derivative factor.
    xi_R=(1.0-xi_M)*relative_R;
    dxi_dseq=(1.0-xi_M)*dxi_dseq;
end
direction = sign(sigma);
eps_R = direction*p.R_H*xi_R;
deps_dsigma = p.R_H*dxi_dseq;
end
