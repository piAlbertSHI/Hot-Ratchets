function [trial,aux] = strain_history_step_1d( ...
    p,state,target_strain,step_duration_s,Tamb)
%SMA.STRAIN_HISTORY_STEP_1D Advance with an explicit strain history.
%
% TARGET_STRAIN contains the endpoint target for each accepted substep and
% STEP_DURATION_S contains its duration.  The quantity is the solver's
% logarithmic strain state (state.eps); callers must document any source
% conversion used to obtain it.  No peak strain, rate, or endpoint is
% inferred from a material parameter or an observation.

if nargin < 5 || isempty(Tamb)
    Tamb = p.Tamb;
end
assert(isvector(target_strain) && ~isempty(target_strain), ...
    'sma:strain_history_step_1d:History', ...
    'An explicit nonempty strain history is required.');
assert(isvector(step_duration_s) && ...
    numel(step_duration_s)==numel(target_strain), ...
    'sma:strain_history_step_1d:History', ...
    'Strain targets and step durations must have equal lengths.');
target_strain = double(target_strain(:));
step_duration_s = double(step_duration_s(:));
assert(all(isfinite(target_strain)) && all(isfinite(step_duration_s)) && ...
    all(step_duration_s>0), ...
    'sma:strain_history_step_1d:History', ...
    'Strain targets must be finite and durations must be positive.');

trial = state;
energy_rate = 0.0;
latent_energy = 0.0;
dissipation_energy = 0.0;
heat_loss_energy = 0.0;
r_phase_energy = 0.0;
strain_work = 0.0;
iterations = 0.0;
residual_norm = 0.0;
thermal_iterations = 0.0;
thermal_temperature_residual = 0.0;
thermal_branch_switches = 0.0;
thermal_cutbacks = 0.0;
monolithic_thermal_solves = 0.0;
thermal_converged = true;

for k = 1:numel(target_strain)
    [trial,local] = sma.coupled_step_1d( ...
        p,trial,1,target_strain(k),Tamb,step_duration_s(k));
    dt = step_duration_s(k);
    energy_rate = energy_rate + get_aux(local,'heat_source_Wkg',0.0)*dt;
    latent_energy = latent_energy + get_aux(local,'latent_specific_Jkg',0.0);
    dissipation_energy = dissipation_energy + ...
        get_aux(local,'dissipation_specific_Jkg',0.0);
    heat_loss_energy = heat_loss_energy + ...
        get_aux(local,'heat_loss_specific_Jkg',0.0);
    r_phase_energy = r_phase_energy + ...
        get_aux(local,'r_phase_latent_specific_Jkg',0.0);
    strain_work = strain_work + get_aux(local,'strain_work_MPa',0.0);
    iterations = iterations + get_aux(local,'iterations',0.0);
    residual_norm = max(residual_norm,get_aux(local,'residual_norm',0.0));
    thermal_iterations = max(thermal_iterations, ...
        get_aux(local,'thermal_iterations',0.0));
    thermal_temperature_residual = max(thermal_temperature_residual, ...
        abs(get_aux(local,'thermomechanical_temperature_residual_K',0.0)));
    thermal_branch_switches = thermal_branch_switches + ...
        get_aux(local,'thermal_branch_switches',0.0);
    thermal_cutbacks = thermal_cutbacks + ...
        get_aux(local,'thermal_cutbacks',0.0);
    monolithic_thermal_solves = monolithic_thermal_solves + ...
        get_aux(local,'monolithic_thermal_solves',0.0);
    thermal_converged = thermal_converged && ...
        logical(get_aux(local,'thermal_converged',true));
    tolerance = max(1.0e-10,10.0*p.numerics.residual_tol* ...
        max(p.Hmax_i,1.0e-4));
    assert(abs(trial.eps-target_strain(k)) <= tolerance, ...
        'sma:strain_history_step_1d:ControlResidual', ...
        'Accepted strain-history substep did not reach its strain target.');
end

aux = local;
aux.control_mode = 'explicit_strain_history';
aux.target_strain = target_strain;
aux.step_duration_s = step_duration_s;
aux.total_duration_s = sum(step_duration_s);
aux.heat_source_Wkg = energy_rate/sum(step_duration_s);
aux.latent_specific_Jkg = latent_energy;
aux.dissipation_specific_Jkg = dissipation_energy;
aux.heat_loss_specific_Jkg = heat_loss_energy;
aux.r_phase_latent_specific_Jkg = r_phase_energy;
aux.strain_work_MPa = strain_work;
aux.iterations = iterations;
aux.residual_norm = residual_norm;
aux.thermal_iterations = thermal_iterations;
aux.thermomechanical_temperature_residual_K = thermal_temperature_residual;
aux.thermal_branch_switches = thermal_branch_switches;
aux.thermal_cutbacks = thermal_cutbacks;
aux.monolithic_thermal_solves = monolithic_thermal_solves;
aux.thermal_converged = thermal_converged;
aux.endpoint_control_residual = trial.eps-target_strain(end);
end

function value = get_aux(aux,name,default_value)
if isfield(aux,name)
    value = aux.(name);
else
    value = default_value;
end
end
