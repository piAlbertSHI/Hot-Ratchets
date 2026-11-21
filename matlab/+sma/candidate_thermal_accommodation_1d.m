function [strain,audit] = candidate_thermal_accommodation_1d( ...
    old_strain,sigma_MPa,rt_activity_weight,T_old_K,T_new_K, ...
    thermal_slope_R_MPaK,Y_R_MPa,moduli_MPa,resistances_MPa)
%SMA.CANDIDATE_THERMAL_ACCOMMODATION_1D Cooling-unlocked fixed-capacity GSM.
% The nonnegative thermodynamic path is available only when forward RT is
% active and the RT thermal force increases during cooling:
%   dTheta = A_R*<s_R*(T_new-T_old)/Y_R>_+,  0<=A_R<=1.
% A_R is supplied by the caller from the RT KKT driving force.  A continuous
% KKT weight avoids both a Heaviside jump at activation and treating the
% positive realmin floor of a log-rate coordinate as physical RT activity.
% The internal strains then follow
%   da_j/dTheta = (sigma-H_j*a_j)/R_j.
% No protocol stage, material name, observation, or endpoint is available.

old_strain=double(old_strain(:));
moduli_MPa=double(moduli_MPa(:));
resistances_MPa=double(resistances_MPa(:));
values=[sigma_MPa,T_old_K,T_new_K,thermal_slope_R_MPaK,Y_R_MPa];
assert(~isempty(old_strain) && numel(old_strain)==numel(moduli_MPa) && ...
    numel(old_strain)==numel(resistances_MPa) && all(isfinite(values)) && ...
    isscalar(rt_activity_weight) && isfinite(rt_activity_weight) && ...
    rt_activity_weight>=0 && rt_activity_weight<=1 && ...
    Y_R_MPa>0, ...
    'sma:candidate_thermal_accommodation_1d:Input', ...
    'Finite thermodynamic data, positive resistance, and matching state widths are required.');

thermal_force_increment=thermal_slope_R_MPaK*(T_new_K-T_old_K)/Y_R_MPa;
unlock_path=rt_activity_weight*max(0,thermal_force_increment);
[strain,base]=sma.candidate_transformation_path_accommodation_1d( ...
    old_strain,sigma_MPa,unlock_path,moduli_MPa,resistances_MPa);

audit=base;
audit.unlock_path=unlock_path;
audit.thermal_force_increment=thermal_force_increment;
audit.rt_activity_weight=rt_activity_weight;
audit.forward_rt_active=rt_activity_weight>0;
audit.temperature_increment_K=T_new_K-T_old_K;
audit.production_connected=false;
end
