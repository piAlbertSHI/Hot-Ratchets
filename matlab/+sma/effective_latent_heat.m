function latent = effective_latent_heat(p,state,trial)
%SMA.EFFECTIVE_LATENT_HEAT Activity- and loading-level-dependent latent heat.

if isfield(trial,'eps_tr_r')
    peak_transformation_strain = max( ...
        abs(state.eps_tr_r),abs(trial.eps_tr_r));
elseif isfield(trial,'h_tr_r')
    peak_transformation_strain = sqrt(2.0/3.0)*max( ...
        norm(state.h_tr_r,'fro'),norm(trial.h_tr_r,'fro'));
else
    peak_transformation_strain = 0.0;
end
load_ratio = min(max(peak_transformation_strain/ ...
    max(p.Hmax_i,eps),0.0),1.0);
activity = max(0.5*(state.zeta_d+trial.zeta_d)-p.zeta_latent0,0.0);
weighted_activity = activity*load_ratio^p.latent_load_power;
load_factor = 1.0-p.latent_load_reduction* ...
    load_ratio^p.latent_load_power;
cycle_factor = 1.0-p.latent_cycle_reduction* ...
    (1.0-exp(-p.latent_cycle_rate*activity))* ...
    load_ratio^p.latent_cycle_power;
latent = load_factor*cycle_factor*(p.latent_f+(p.latent-p.latent_f)* ...
    exp(-p.lambda_latent*weighted_activity));
end
