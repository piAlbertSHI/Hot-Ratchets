function authority = protocol_authority(case_name)
%SMA.PROTOCOL_AUTHORITY Source facts and evidence gates for literature cards.
%
% This metadata is deliberately separate from the numerical surrogate path
% returned by SMA.LITERATURE_PROTOCOL.  It records what the papers actually
% specify, which control history is still missing, and whether a result may be
% used as experimental validation.  It does not insert observations into a
% constitutive update and it does not invent an unreported force history.

if nargin < 1 || isempty(case_name)
    error('sma:protocol_authority:Case', ...
        'A literature case name is required.');
end
name = lower(char(case_name));
case_id = sma.case_id(name);

authority = struct( ...
    'case_name',name, ...
    'case_id',case_id, ...
    'source_family','', ...
    'experiment_family','', ...
    'loading_control','', ...
    'hold_control','', ...
    'unloading_control','', ...
    'recovery_control','', ...
    'source_loading_rate_1_per_s',NaN, ...
    'source_rate_owner','', ...
    'protocol_card_id','', ...
    'source_hold_after_load_s',NaN, ...
    'source_hold_after_unload_s',NaN, ...
    'source_peak_stress_MPa',NaN, ...
    'source_test_temperature_K',NaN, ...
    'source_total_unloading_time_s',NaN, ...
    'source_strain_measure','', ...
    'source_stress_measure','', ...
    'source_area_basis','', ...
    'source_area_measure_identified',false, ...
    'source_measure_mapping_status','', ...
    'execution_strain_measure','', ...
    'execution_stress_measure','', ...
    'execution_area_basis','', ...
    'source_force_history_complete',false, ...
    'source_strain_history_complete',false, ...
    'surrogate_control_mode','', ...
    'surrogate_stress_rate_MPa_per_s',NaN, ...
    'surrogate_mapping_status','', ...
    'execution_control_mode','', ...
    'execution_rate_1_per_s',NaN, ...
    'execution_is_source_authoritative',false, ...
    'experimental_acceptance_allowed',false);

switch case_id
    case 1
        authority.source_family = 'Zhang2022_NiTi';
        authority.experiment_family = 'cyclic_tension_unloading_eCE';
        authority.protocol_card_id = ['zhang2022_niti_' name];
        authority.loading_control = 'displacement_or_current_length_strain';
        authority.hold_control = 'displacement_at_peak';
        authority.unloading_control = 'force_to_zero_stress';
        authority.recovery_control = 'zero_stress_hold';
        authority.source_loading_rate_1_per_s = 1.0e-3;
        authority.source_rate_owner = 'Zhang_cyclic_loading_history';
        authority.source_total_unloading_time_s = ...
            zhang_unloading_time_from_name(name);
        authority.source_strain_measure = 'current_length_normalized_(l-l0)/l';
        authority.source_stress_measure = 'reported_stress_area_operator_unresolved';
        authority.source_area_basis = 'unreported_reference_or_current';
        authority.source_area_measure_identified = false;
        authority.source_measure_mapping_status = ...
            'source_current_length_and_stress_area_not_closed';
        authority.execution_strain_measure = 'additive_small_strain_work_conjugate';
        authority.execution_stress_measure = 'small_strain_work_conjugate_stress';
        authority.execution_area_basis = 'model_internal_reference_measure';
        % The paper gives terminal zero stress and total duration, but no
        % force setpoint trajectory/rate.  The current runner's 80 MPa/s path
        % is therefore a declared surrogate, not a source-authoritative card.
        authority.source_force_history_complete = false;
        authority.source_strain_history_complete = true;
        authority.surrogate_control_mode = 'stress_rate_surrogate';
        authority.surrogate_stress_rate_MPa_per_s = 80.0;
        authority.surrogate_mapping_status = ...
            'conditional_surrogate_not_source_authoritative';
        authority.execution_control_mode = 'mixed_strain_loading_stress_surrogate_unload';
        authority.execution_rate_1_per_s = NaN;
        authority.execution_is_source_authoritative = false;
    case {2,3,4}
        authority.source_family = 'Yuan2020_CuAlMn';
        authority.experiment_family = 'direct_eCE_load_unload';
        if case_id==2
            authority.protocol_card_id = 'yuan2020_direct_ece_annealed_130um';
            authority.source_peak_stress_MPa = 150.0;
        elseif case_id==3
            authority.protocol_card_id = 'yuan2020_direct_ece_drawn_130um';
            authority.source_peak_stress_MPa = 300.0;
        else
            authority.protocol_card_id = 'yuan2020_direct_ece_drawn_80um';
            authority.source_peak_stress_MPa = 300.0;
        end
        authority.loading_control = 'strain_rate_controlled';
        authority.hold_control = 'stress_or_load_hold';
        authority.unloading_control = 'strain_rate_controlled';
        authority.recovery_control = 'zero_stress_hold';
        authority.source_loading_rate_1_per_s = 0.15;
        authority.source_rate_owner = 'Yuan2020_direct_eCE_load_unload';
        authority.source_hold_after_load_s = 15.0;
        authority.source_hold_after_unload_s = 15.0;
        authority.source_test_temperature_K = 303.0;
        authority.source_total_unloading_time_s = NaN;
        authority.source_strain_measure = 'engineering_(l-l0)/l0';
        authority.source_stress_measure = 'reported_stress_area_operator_unresolved';
        authority.source_area_basis = 'unreported_reference_or_current';
        authority.source_area_measure_identified = false;
        authority.source_measure_mapping_status = ...
            'source_engineering_strain_declared_stress_area_not_closed';
        authority.execution_strain_measure = 'logarithmic_solver_state_with_engineering_input';
        authority.execution_stress_measure = 'small_strain_work_conjugate_stress';
        authority.execution_area_basis = 'model_internal_reference_measure';
        % Fig. 5/6 state the strain rate and 15 s holds, but the complete
        % peak-strain history is not tabulated for every specimen.  The same
        % paper's superelastic force-controlled tests are a different card.
        authority.source_force_history_complete = false;
        authority.source_strain_history_complete = false;
        authority.surrogate_control_mode = 'stress_control_surrogate';
        authority.surrogate_stress_rate_MPa_per_s = NaN;
        authority.surrogate_mapping_status = ...
            'conditional_surrogate_not_source_authoritative';
        authority.execution_control_mode = 'stress_control_surrogate';
        % The numerical protocol below intentionally does not claim to
        % execute the source rate until a strain-time controller trace is
        % available.  Keep this NaN rather than silently equating a 1 s
        % stress ramp with Yuan's 0.15 s^-1 strain rate.
        authority.execution_rate_1_per_s = NaN;
        authority.execution_is_source_authoritative = false;
    case {5,7}
        authority.source_family = 'ImranZhang2020_secondary_review';
        authority.protocol_card_id = [authority.source_family '_' name];
        authority.experiment_family = 'secondary_source_trend_benchmark';
        authority.loading_control = 'not_source_complete';
        authority.hold_control = 'not_source_complete';
        authority.unloading_control = 'not_source_complete';
        authority.recovery_control = 'not_source_complete';
        authority.surrogate_mapping_status = ...
            'secondary_trend_only_not_experimental_validation';
        authority.execution_control_mode = 'secondary_trend_surrogate';
        authority.source_strain_measure = 'secondary_source_measure_unresolved';
        authority.source_stress_measure = 'secondary_source_measure_unresolved';
        authority.source_area_basis = 'unresolved';
        authority.source_area_measure_identified = false;
        authority.source_measure_mapping_status = 'not_source_authoritative';
        authority.execution_strain_measure = 'surrogate_only';
        authority.execution_stress_measure = 'surrogate_only';
        authority.execution_area_basis = 'surrogate_only';
        authority.execution_rate_1_per_s = NaN;
        authority.execution_is_source_authoritative = false;
    otherwise
        error('sma:protocol_authority:UnsupportedCase', ...
            'No source-authority metadata for case ID %d.',case_id);
end
end

function value = zhang_unloading_time_from_name(name)
tokens = regexp(name,'niti_s(4|6|8|10|12)$','tokens','once');
if isempty(tokens)
    value = NaN;
else
    value = str2double(tokens{1});
end
end
