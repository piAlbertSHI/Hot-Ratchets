function value = r_phase_latent_specific( ...
    p,xi_old,seq_old,T_old,training_old, ...
    xi_new,seq_new,T_new,training_new,xi_R_old,xi_R_new)
%SMA.R_PHASE_LATENT_SPECIFIC Reduced A/R latent-heat increment.

value = 0.0;
if ~p.enable_r_phase || p.R_latent_scale<=0.0
    return;
end

if p.enable_explicit_r_phase_state
    assert(nargin>=11,'sma:r_phase_latent_specific:CommittedState', ...
        'Explicit RT latent heat requires committed R fractions.');
    fraction_old=xi_R_old;
    fraction_new=xi_R_new;
else
    fraction_old = (1.0-min(max(xi_old,0.0),1.0))* ...
        sma.r_phase_fraction(p,seq_old,T_old,training_old);
    fraction_new = (1.0-min(max(xi_new,0.0),1.0))* ...
        sma.r_phase_fraction(p,seq_new,T_new,training_new);
end
Tmean = 0.5*(T_old+T_new);
latent_Jkg = -Tmean*p.R_entropy*1.0e6/p.rho;
value = p.R_latent_scale*latent_Jkg*(fraction_new-fraction_old);
end
