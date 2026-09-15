function [h_R,xi_R] = r_phase_strain_3d( ...
    p,tau,T,training_activity,xi_M,xi_R_old,xi_M_old)
%SMA.R_PHASE_STRAIN_3D Reversible J2 R-phase transformation strain.

if nargin<4
    training_activity = 0.0;
end
if nargin<5
    xi_M=0.0;
end
[seq,dev] = sma.equivalent_stress(tau);
if p.enable_r_phase_return_map
    assert(nargin>=7,'sma:r_phase_strain_3d:CommittedState', ...
        'The RT return map requires committed R and M fractions.');
    [~,xi_R] = sma.r_phase_return_map( ...
        p,xi_M,seq,T,training_activity,xi_M_old,xi_R_old);
elseif p.enable_explicit_r_phase_state
    [~,xi_R] = sma.explicit_phase_fractions( ...
        p,xi_M,seq,T,training_activity);
else
    % Return the absolute R inventory on every branch.  A relative R fraction
    % is not conjugate to the mixture strain and was inconsistent with the
    % latent-heat routine's (1-xi_M) weighting.
    xi_R = (1.0-xi_M) * sma.r_phase_fraction(p,seq,T,training_activity);
end
if seq>p.numerics.direction_tol
    N = 1.5*dev/seq;
else
    N = zeros(3);
end
h_R = p.R_H*xi_R*N;
end
