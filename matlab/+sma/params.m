function p = params(name)
%SMA.PARAMS Unique material parameter definitions for the standard model.
% Stress MPa; temperature K; strain dimensionless; energy J/kg.

p = base_parameters();
switch lower(char(name))
    case {'niti_zhang2022', 'niti'}
        p.name = 'NiTi Zhang 2022 S8';
        p.case_id = 1;
        p.T0 = 301.0;
        p.Tamb = 301.0;
        % Direct anchors from Zhang et al. (2022), Table 2 and Figure 1.
        p.EA = 48.4e3;
        p.EM = 25.4e3;
        p.nuA = 0.30;
        p.nuM = 0.30;
        p.alphaA = 22e-6;
        p.alphaM = 22e-6;
        p.Ms = 205.0;
        % Mf/As/Af close the reduced MT branch after the reversible
        % algebraic R-phase contribution has been separated.
        p.Mf = 185.0;
        p.As = 222.0;
        p.Af = 255.0;
        p.CM = 6.0;
        p.CA = 5.0;
        p.rho_delta_s = -0.23;
        p.Hmin = 0.0035;
        p.Hmax_i = 0.0915;
        p.Hmax_f = 0.0910;
        p.kt = 0.0330;
        p.sigma_crit = 80.0;
        p.lambda_H = 0.005;
        p.Cp1 = 0.0054;
        p.Cp2 = 1.20;
        p.sigma_b = 0.0;
        p.lambda_beta = 0.18;
        % Effective one-fraction heat closure. The paper's DSC heat
        % cannot be copied directly because RT and MT heat are collapsed.
        p.latent = 6.20e3;
        p.latent_f = 4.60e3;
        p.lambda_latent = 0.20;
        p.latent_load_power = 20.0;
        p.latent_load_reduction = 0.16;
        p.latent_cycle_power = 8.0;
        p.latent_cycle_reduction = 0.50;
        p.latent_cycle_rate = 0.20;
        p.zeta_latent0 = 2.0;
        p.heat_loss_rate = 12.0;
        p.hold_heat_multiplier = 5.0;
        p.unload_substeps = 32;
        p.unload_heat_multiplier = 0.10;
        p.unload_max_cutbacks = 16;
        p.reference_stress = 500.0;
        % The paper reports Yf0=21.7 MPa. The one-fraction closure folds
        % the omitted R-phase resistance into the effective values below.
        p.Yf = 36.0;
        p.Yf_f = 34.0;
        p.Yr = 1.1;
        p.Yr_f = 0.5;
        p.lambda_Y = 0.50;
        p.zeta_Y0 = 0.0;
        % Zhang et al. (2022), Table 2. The R phase is reversible and
        % algebraic in this reduced model, so it adds no history variable.
        p.enable_r_phase = true;
        p.R_H = 0.0035;
        p.R_entropy = -0.06;
        p.R_tau0 = 1.3;
        p.R_dissipation = 0.1;
        p.R_hardening = 1.4;
        p.R_Tref = 301.0;
        p.R_latent_scale = 0.20;
        p.R_training_gain = 4.0;
        p.R_training_rate = 0.30;
        % Zhang et al. (2022), Table 2: A/R-phase plasticity (ARP).
        p.enable_arp = true;
        p.arp_tau0 = 600.0;
        p.arp_tau0_f = 150.0;
        p.arp_tau_fraction_start = 0.25;
        p.arp_tau_fraction_width = 0.35;
        p.arp_H = 150.0e3;
        p.arp_matrix_power = 1.0;
        % The linear term is absorbed by the endpoint closure. The reported
        % ninth-order forward and reverse terms are retained explicitly.
        p.enable_mt_nonlinear_hardening = true;
        p.mt_H = 0.0;
        p.mt_hfor = 72.0;
        p.mt_hfor_training_extra = 26.0;
        p.mt_hfor_training_scale = 10.0;
        p.mt_hfor_training_power = 1.0;
        p.mt_hfor_transition_dip = 2.4;
        p.mt_hfor_transition_width = 0.14;
        p.mt_hfor_shoulder = 1.3;
        p.mt_hfor_shoulder_start = 0.29;
        p.mt_hfor_shoulder_width = 0.06;
        p.mt_hfor_shoulder_memory_scale = 0.030;
        p.mt_hfor_shoulder_memory_power = 8.0;
        p.mt_hfor_shoulder_training_scale = 10.0;
        p.mt_hfor_shoulder_training_power = 1.0;
        % High-cycle/high-transformation hardening closes the S12 cycle-20
        % loading endpoint without shifting lower-strain or first-cycle loops.
        p.mt_hfor_activity = 7000.0;
        p.mt_hfor_activity_power = 20.0;
        p.mt_hrev = 20.0;
        p.mt_hrev_span = 8.0;
        p.mt_hrev_span_load_power = 4.0;
        p.mt_hrev_finish = 8.0;
        p.mt_hrev_finish_low = 12.0;
        p.mt_hrev_finish_low_load_power = 2.0;
        p.mt_hrev_finish_low_power = 2.0;
        p.mt_hrev_finish_low_training_scale = 1.0;
        p.mt_hrev_finish_low_training_power = 1.0;
        p.mt_hrev_finish_load_power = 2.0;
        p.mt_hrev_finish_power = 2.0;
        p.mt_reverse_memory_shift = 40.0;
        p.mt_reverse_memory_scale = 0.040;
        p.mt_reverse_memory_power = 8.0;
        p.mt_n = 9.0;
        % Zhang et al. (2022), Eqs. (54)-(57): Y=Yf-Yp,
        % Yf=Yf0+d*xi_recent^n and dYp=aY*(cY-Yp)*|dxi|. The reduced
        % cycle-frozen activity has its own calibrated time scale.
        p.enable_mt_dissipation_evolution = true;
        p.mt_Y_d = 2.8;
        p.mt_Y_n = 5.0;
        % Source-card literals for Zhang (2022) Eq. (39).  These are kept
        % separate from the reduced one-fraction Yf closure above.
        p.zhang_Yf0_MPa = 21.7;
        p.zhang_Y_d_MPa = 2.8;
        p.zhang_Y_n_xi = 5.0;
        p.zhang_m_tr = 200.0;
        p.zhang_reference_rate_per_s = 1.0;
        p.mt_Yp_sat = 13.0;
        p.mt_Yp_high_sat = 45.0;
        p.mt_Yp_memory_scale = 0.040;
        p.mt_Yp_memory_power = 8.0;
        p.mt_Yp_activity_scale = 35.0;
        p.mt_Yp_activity_power = 6.0;
        p.mt_Yp_rate = 0.02;
        p.mt_Yp_reverse_scale = 10.2;
        p.mt_Yp_reverse_training_decay = 0.12;
        p.mt_Yp_reverse_memory_scale = 0.015;
        p.mt_Yp_reverse_memory_power = 8.0;
        p.mt_Yp_use_cumulative_activity = false;
        p.mt_training_reference_fraction = 0.25;
        p.mt_training_cycle_increment = 2.0;
        p.mt_training_completion_tol = 3.0e-4;
        % Reduced closure of Zhang et al. (2022), Eqs. (53b) and (58).
        % The saturation envelope is calibrated once from Fig. 14.
        p.enable_residual_martensite = true;
        p.residual_xi_a = 0.60;
        p.residual_xi_offset = 0.006;
        p.residual_xi_peak_gain = 0.034;
    case {'cualmn_annealed', 'annealed'}
        p = cu_base();
        p.name = 'Cu-Al-Mn annealed 130 um';
        p.case_id = 2;
        p.EA = 38e3;
        p.EM = 24e3;
        p.Ms = 250.0;
        p.Mf = 236.0;
        p.As = 274.0;
        p.Af = 291.0;
        p.CM = 1.05;
        p.CA = 1.05;
        p.rho_delta_s = -0.0672;
        p.Hmax_i = 0.064;
        p.Hmax_f = 0.060;
        p.lambda_H = 0.010;
        p.zeta_H0 = 24.0;
        p.Cp1 = 0.0;
        p.Cp2 = 0.90;
        p.rho = 7500.0;
        p.cp = 297.0;
        p.latent = 1.515e3;
        p.latent_f = 0.900e3;
        p.lambda_latent = 0.080;
        p.zeta_latent0 = 24.0;
        p.reference_stress = 150.0;
    case {'cualmn_drawn130', 'drawn130'}
        p = cu_base();
        p.name = 'Cu-Al-Mn cold-drawn 130 um';
        p.case_id = 3;
        p.EA = 38e3;
        p.EM = 25e3;
        p.Ms = 242.0;
        p.Mf = 229.0;
        p.As = 264.0;
        p.Af = 279.0;
        p.CM = 2.25;
        p.CA = 2.25;
        p.rho_delta_s = -0.1305;
        p.Hmax_i = 0.058;
        p.Hmax_f = 0.055;
        p.lambda_H = 0.002;
        % The Fig. 6(b) envelope stays nearly flat for roughly 25 cycles.
        % Use a continuous defect-incubation threshold rather than a
        % cycle-number switch; degradation starts only after sufficient
        % transformation activity has accumulated.
        p.zeta_H0 = 56.0;
        p.Cp1 = 0.010;
        p.Cp2 = 0.70;
        p.rho = 7500.0;
        p.cp = 410.0;
        p.latent = 2.10e3;
        p.latent_f = 1.45e3;
        p.lambda_latent = 0.030;
        p.zeta_latent0 = 56.0;
        p.reference_stress = 300.0;
    case {'cualmn_drawn80', 'drawn80'}
        p = cu_base();
        p.name = 'Cu-Al-Mn cold-drawn 80 um';
        p.case_id = 4;
        p.EA = 40e3;
        p.EM = 27e3;
        p.Ms = 242.0;
        p.Mf = 229.0;
        p.As = 264.0;
        p.Af = 279.0;
        p.CM = 2.25;
        p.CA = 2.25;
        p.rho_delta_s = -0.1305;
        p.Hmax_i = 0.040;
        p.Hmax_f = 0.0395;
        p.lambda_H = 0.001;
        p.Cp1 = 0.002;
        p.Cp2 = 0.50;
        p.rho = 7500.0;
        p.cp = 410.0;
        p.latent = 1.50e3;
        p.latent_f = 1.46e3;
        p.lambda_latent = 0.010;
        p.reference_stress = 300.0;
    case {'nimnga_textured_review', 'nimnga', 'review'}
        p = nimnga_base();
        p.name = 'Textured Ni-Mn-Ga review Figure 5 benchmark';
        p.case_id = 5;
    case {'nimnga_review_fig4', 'review_fig4'}
        p = nimnga_base();
        p.name = 'Textured Ni-Mn-Ga review Figure 4 benchmark';
        p.case_id = 7;
        p.T0 = 300.0;
        p.Tamb = 300.0;
        p.Ms = 285.0;
        p.Mf = 275.0;
        p.As = 287.0;
        p.Af = 297.0;
        p.latent = 2.05e3;
        p.latent_f = 2.05e3;
        p.Hmax_f = 0.040;
        p.lambda_H = 0.0;
        p.Cp1 = 0.0;
        p.reference_stress = 65.0;
    case {'qian2019_se21', 'se21'}
        p.name = 'Qian 2019 SE21 NiTi';
        p.case_id = 6;
        p.T0 = 303.15;
        p.Tamb = 303.15;
        p.EA = 30e3;
        p.EM = 30e3;
        p.Ms = 258.15;
        p.Mf = 248.15;
        p.As = 268.15;
        p.Af = 278.15;
        p.CM = 7.3;
        p.CA = 7.6;
        p.cp = 500.0;
        p.alphaA = 0.0;
        p.alphaM = 0.0;
        p.Hmin = 0.0;
        p.Hmax_i = 0.030;
        p.Hmax_f = 0.030;
        p.kt = 0.030;
        p.sigma_crit = 0.0;
        p.lambda_H = 0.0;
        p.Cp1 = 0.0;
        p.Cp2 = 0.5;
        p.sigma_b = 0.0;
        p.lambda_beta = 0.1;
        p.latent = 10.0e3;
        p.heat_loss_rate = 0.0;
        p.reference_stress = 450.0;
    otherwise
        error('sma:params:UnknownSet', 'Unknown parameter set: %s', char(name));
end
p = validate_parameters(p);
end

function p = base_parameters()
p = struct();
p.name = 'base';
p.case_id = 0;
p.T0 = 303.15;
p.Tamb = 303.15;
p.EA = 45e3;
p.EM = 28e3;
p.nuA = 0.33;
p.nuM = 0.33;
p.alphaA = 1.0e-5;
p.alphaM = 1.0e-5;
p.Ms = 285.0;
p.Mf = 265.0;
p.As = 280.0;
p.Af = 300.0;
p.CM = 6.0;
p.CA = 6.0;
p.Hmin = 0.0;
p.Hmax_i = 0.040;
p.Hmax_f = 0.040;
p.kt = 0.020;
p.sigma_crit = 0.0;
p.lambda_H = 0.0;
p.zeta_H0 = 0.0;
p.Cp1 = 0.005;
p.Cp2 = 0.50;
p.sigma_b = 20.0;
p.lambda_beta = 0.10;
p.rho = 6500.0;
p.cp = 410.0;
p.latent = 2.1e3;
p.latent_f = NaN;
p.lambda_latent = 0.0;
p.latent_load_power = 1.0;
p.latent_load_reduction = 0.0;
p.latent_cycle_power = 1.0;
p.latent_cycle_reduction = 0.0;
p.latent_cycle_rate = 0.0;
p.zeta_latent0 = 0.0;
p.taylor_quinney = 0.10;
p.heat_loss_rate = 8.0;
p.hold_heat_multiplier = 1.0;
p.hold_target_tol = 1.0e-9;
p.unload_substeps = 1;
p.unload_heat_multiplier = 1.0;
p.unload_max_cutbacks = 4;
p.reference_stress = 200.0;
p.rho_delta_s = NaN;
p.rho_delta_c = 0.0;
p.rho_delta_u = 0.0;
p.Yf = 0.0;
p.Yr = 0.0;
p.Yf_f = NaN;
p.Yr_f = NaN;
p.lambda_Y = 0.0;
p.zeta_Y0 = 0.0;
% Continuous two-population resistance relaxation. NaN rates are resolved
% to lambda_Y after each material card is assembled, which reproduces the
% legacy single exponential exactly without a mechanism switch.
p.Y_activity_fast_weight = 1.0;
p.Y_activity_fast_rate = NaN;
p.Y_activity_slow_rate = NaN;
p.Y_activity_partial_power = 0.0;
% Optional reduced active-population moment.  When enabled, the same
% resistance law is driven by the committed traversal activity zeta_p
% instead of the total transformed-volume path zeta_d.
p.Y_activity_use_population_path = false;
% Versioned explicit threshold-population state.  The production mixing
% coefficient remains zero until the v2 return-map and cross-material gates
% pass; bin count and width are numerical resolution parameters.
p.Y_activity_population_mix = 0.0;
p.Y_activity_population_bins = 8;
p.Y_activity_population_width = 1.0/8.0;
% Direct Zhang (2022) Eqs. (48)-(60) history.  A positive coupling uses
% only the directly closed Yp equation in the transformation surface; the
% remaining microstructural variables are committed diagnostics.
p.zhang_direct_Yp_coupling = 0.0;
p.zhang_direct_Omega_coupling = 0.0;
p.enable_explicit_r_phase_state = false;
p.enable_r_phase_return_map = false;
p.enable_surface_driven_branch_selection = false;
% Zhang Eqs. (6) and (37-b): use the current explicit A/R inventory in the
% MT transformation-strain magnitude and structural-entropy difference.
% The shadow switch is neutral until the two-surface gate is passed.
p.enable_composition_weighted_mt = false;
% Zhang Eqs. (61)-(63) source-resistance shadow.  Neutral by default so the
% calibrated reduced card remains the production path until a real gate passes.
p.enable_zhang_source_mt_resistance = false;
p.zhang_tau0_AR_M = 50.0;
p.zhang_H_AR_M = 10.0;
p.zhang_hfor_AR_M = 46.0;
p.zhang_hrev_AR_M = 20.0;
p.zhang_ntr_AR_M = 9.0;
% Optional converged thermo-mechanical fixed point.  The default preserves
% the existing two-pass staggered update until the cross-amplitude gate is
% passed.  Every fixed-point probe is evaluated from the committed state.
p.enable_iterated_thermomechanical_coupling = false;
p.thermomechanical_coupling_max_iterations = 40;
p.thermomechanical_coupling_temperature_tol_K = 1.0e-7;
p.thermomechanical_coupling_relaxation = 1.0;
p.thermomechanical_coupling_max_cutbacks = 8;
p.enable_r_phase = false;
p.R_H = 0.0;
p.R_entropy = 0.0;
p.R_tau0 = 0.0;
p.R_dissipation = 0.0;
p.R_hardening = 1.0;
p.R_Tref = NaN;
p.R_max_fraction = 1.0;
p.R_latent_scale = 0.0;
p.R_training_gain = 0.0;
p.R_training_rate = 0.0;
p.enable_arp = false;
p.arp_tau0 = 0.0;
p.arp_tau0_f = 0.0;
p.arp_tau_fraction_start = 1.0;
p.arp_tau_fraction_width = 1.0;
p.arp_H = 1.0;
p.arp_matrix_power = 1.0;
p.enable_mt_nonlinear_hardening = false;
p.mt_H = 0.0;
p.mt_hfor = 0.0;
p.mt_hfor_training_extra = 0.0;
p.mt_hfor_training_scale = 1.0;
p.mt_hfor_training_power = 1.0;
p.mt_hfor_transition_dip = 0.0;
p.mt_hfor_transition_width = 1.0;
p.mt_hfor_shoulder = 0.0;
p.mt_hfor_shoulder_start = 1.0;
p.mt_hfor_shoulder_width = 1.0;
p.mt_hfor_shoulder_memory_scale = 1.0;
p.mt_hfor_shoulder_memory_power = 1.0;
p.mt_hfor_shoulder_training_scale = 1.0;
p.mt_hfor_shoulder_training_power = 1.0;
p.mt_hfor_activity = 0.0;
p.mt_hfor_activity_power = 1.0;
p.mt_hrev = 0.0;
p.mt_hrev_span = 0.0;
p.mt_hrev_span_load_power = 1.0;
p.mt_hrev_finish = 0.0;
p.mt_hrev_finish_low = 0.0;
p.mt_hrev_finish_low_load_power = 1.0;
p.mt_hrev_finish_low_power = 1.0;
p.mt_hrev_finish_low_training_scale = 1.0;
p.mt_hrev_finish_low_training_power = 1.0;
p.mt_hrev_finish_load_power = 0.0;
p.mt_hrev_finish_power = 1.0;
p.mt_reverse_memory_shift = 0.0;
p.mt_reverse_memory_scale = 1.0;
p.mt_reverse_memory_power = 1.0;
p.mt_n = 1.0;
p.enable_mt_dissipation_evolution = false;
p.mt_Y_d = 0.0;
p.mt_Y_n = 1.0;
p.mt_Yp_sat = 0.0;
p.mt_Yp_high_sat = 0.0;
p.mt_Yp_memory_scale = 1.0;
p.mt_Yp_memory_power = 1.0;
p.mt_Yp_activity_scale = 1.0;
p.mt_Yp_activity_power = 1.0;
p.mt_Yp_rate = 0.0;
p.mt_Yp_reverse_scale = 1.0;
p.mt_Yp_reverse_training_decay = 0.0;
p.mt_Yp_reverse_memory_scale = 1.0;
p.mt_Yp_reverse_memory_power = 1.0;
p.mt_Yp_use_cumulative_activity = false;
p.mt_training_reference_fraction = 1.0;
p.mt_training_cycle_increment = 0.0;
p.mt_training_completion_tol = 1.0e-8;
p.enable_residual_martensite = false;
p.residual_xi_a = 0.0;
p.residual_xi_offset = 0.0;
p.residual_xi_peak_gain = 0.0;
p.min_temperature = 200.0;
p.max_temperature = 500.0;
p.default_dt = 0.20;
p.enable_rate_regularization = false;
p.phase_loading_rate_1_per_s = inf;
p.phase_unloading_rate_1_per_s = inf;
p.phase_overdrive_width_MPa = 1.0;
p.phase_overdrive_activation_enabled = false;
p.phase_overdrive_activation_shape = 'cubic_smoothstep';
p.phase_rate_multiplier_mode = ...
    'none_qian_eq3_direction_coefficient_exact';
p.phase_kinetics_model = 'rate_independent';
p.numerics = struct( ...
    'surface_tol', 1.0e-8, ...
    'residual_tol', 2.0e-6, ...
    'step_tol', 1.0e-9, ...
    'max_iterations', 40, ...
    'max_line_search', 12, ...
    'jacobian_step', 2.0e-6, ...
    'direction_tol', 1.0e-10, ...
    'xi_tol', 1.0e-10, ...
    'stress_solve_bound', 2500.0, ...
    'thermal_iterations', 3, ...
    'cutback_on_failure', true);
end

function p = cu_base()
p = base_parameters();
p.T0 = 303.0;
p.Tamb = 303.0;
p.Hmin = 0.001;
p.kt = 0.018;
p.sigma_crit = 20.0;
p.sigma_b = 0.0;
p.lambda_beta = 0.08;
p.heat_loss_rate = 82.0;
p.default_dt = 0.25;
end

function p = nimnga_base()
p = base_parameters();
p.T0 = 327.0;
p.Tamb = 327.0;
p.EA = 45e3;
p.EM = 27e3;
p.Ms = 295.0;
p.Mf = 285.0;
p.As = 320.0;
p.Af = 327.0;
p.CM = 2.4;
p.CA = 2.6;
p.Hmin = 0.038;
p.Hmax_i = 0.040;
p.Hmax_f = 0.0395;
p.kt = 0.040;
p.sigma_crit = 5.0;
p.lambda_H = 0.0;
p.Cp1 = 0.0;
p.Cp2 = 0.25;
p.sigma_b = 1.0;
p.lambda_beta = 0.05;
p.latent = 2.15e3;
p.latent_f = 2.15e3;
p.unload_substeps = 8;
p.unload_heat_multiplier = 0.10;
p.heat_loss_rate = 82.0;
p.reference_stress = 100.0;
end

function p = validate_parameters(p)
assert(p.EA > 0 && p.EM > 0, 'Elastic moduli must be positive.');
assert(p.Mf < p.Ms, 'Mf must be below Ms.');
assert(p.As < p.Af, 'As must be below Af.');
assert(p.Hmax_i >= p.Hmin && p.Hmax_f >= p.Hmin, ...
    'Hmax must not be below Hmin.');
assert(p.cp > 0 && p.rho > 0, 'Density and heat capacity must be positive.');
if ~isfinite(p.rho_delta_s)
    % MPa/K. The entropy jump is negative for A->M. The phase-endpoint
    % hardening functions then recover Ms/Mf/As/Af at zero stress.
    p.rho_delta_s = -p.Hmax_i * 0.5 * (p.CM + p.CA);
end
if ~isfinite(p.latent_f)
    p.latent_f = p.latent;
end
if ~isfinite(p.Yf_f)
    p.Yf_f = p.Yf;
end
if ~isfinite(p.Yr_f)
    p.Yr_f = p.Yr;
end
if ~isfinite(p.R_Tref)
    p.R_Tref = p.T0;
end
if ~isfinite(p.Y_activity_fast_rate)
    p.Y_activity_fast_rate = p.lambda_Y;
end
if ~isfinite(p.Y_activity_slow_rate)
    p.Y_activity_slow_rate = p.lambda_Y;
end
assert(p.latent_f>0 && p.latent_f<=p.latent, ...
    'Final latent heat must be positive and no larger than its initial value.');
assert(p.lambda_latent>=0 && p.latent_load_power>=0 && ...
    p.latent_load_reduction>=0 && p.latent_load_reduction<1, ...
    'Latent-heat degradation parameters must be nonnegative.');
assert(p.latent_cycle_power>=0 && p.latent_cycle_reduction>=0 && ...
    p.latent_cycle_reduction<1 && p.latent_cycle_rate>=0, ...
    'Cycle-dependent entropy-loss parameters must be admissible.');
assert(p.Y_activity_fast_weight>=0 && p.Y_activity_fast_weight<=1 && ...
    p.Y_activity_fast_rate>=0 && p.Y_activity_slow_rate>=0 && ...
    p.Y_activity_partial_power>=0, ...
    'Resistance activity-spectrum parameters must be admissible.');
assert((islogical(p.Y_activity_use_population_path) || ...
    isnumeric(p.Y_activity_use_population_path)) && ...
    isscalar(p.Y_activity_use_population_path), ...
    'Population-path selector must be a scalar logical value.');
assert(p.Y_activity_population_mix>=0 && ...
    p.Y_activity_population_mix<=1 && ...
    p.Y_activity_population_bins>=2 && ...
    p.Y_activity_population_bins==round(p.Y_activity_population_bins) && ...
    p.Y_activity_population_width>0, ...
    'Threshold-population parameters must be admissible.');
assert(p.zhang_direct_Yp_coupling>=0 && ...
    p.zhang_direct_Yp_coupling<=1, ...
    'Direct Zhang Yp coupling must lie in [0,1].');
assert(p.zhang_direct_Omega_coupling>=0 && ...
    p.zhang_direct_Omega_coupling<=1, ...
    'Direct Zhang Omega coupling must lie in [0,1].');
assert(isscalar(p.enable_explicit_r_phase_state), ...
    'Explicit R-phase state switch must be scalar.');
assert(isscalar(p.enable_r_phase_return_map) && ...
    (~p.enable_r_phase_return_map || p.enable_explicit_r_phase_state), ...
    'The RT return map requires the explicit R-phase state.');
assert(isscalar(p.enable_surface_driven_branch_selection), ...
    'Surface-driven branch selection switch must be scalar.');
assert(isscalar(p.enable_composition_weighted_mt) && ...
    (~p.enable_composition_weighted_mt || ...
    (p.enable_explicit_r_phase_state && p.enable_r_phase_return_map)), ...
    ['Composition-weighted MT requires the explicit conservative R state ' ...
    'and its independent return map.']);
assert(isscalar(p.enable_zhang_source_mt_resistance) && ...
    p.zhang_tau0_AR_M>=0 && p.zhang_H_AR_M>=0 && ...
    p.zhang_hfor_AR_M>=0 && p.zhang_hrev_AR_M>=0 && ...
    p.zhang_ntr_AR_M>=1, ...
    'Zhang source MT resistance controls are not admissible.');
assert(isscalar(p.enable_iterated_thermomechanical_coupling) && ...
    p.thermomechanical_coupling_max_iterations>=2 && ...
    p.thermomechanical_coupling_temperature_tol_K>0 && ...
    p.thermomechanical_coupling_relaxation>0 && ...
    p.thermomechanical_coupling_relaxation<=1 && ...
    p.thermomechanical_coupling_max_cutbacks>=0 && ...
    p.thermomechanical_coupling_max_cutbacks== ...
        round(p.thermomechanical_coupling_max_cutbacks), ...
    'Iterated thermo-mechanical coupling controls are not admissible.');
assert(p.hold_heat_multiplier>=1 && p.hold_target_tol>=0, ...
    'Hold heat-transfer controls must be admissible.');
assert(p.unload_substeps>=1, ...
    'Fast-unloading substep count must be at least one.');
assert(p.unload_heat_multiplier>=0, ...
    'Fast-unloading heat-transfer multiplier must be nonnegative.');
assert(p.unload_max_cutbacks>=0, ...
    'Fast-unloading cutback count must be nonnegative.');
assert(p.Yf>=0 && p.Yr>=0 && p.Yf_f>=0 && p.Yr_f>=0, ...
    'Transformation resistances must be nonnegative.');
assert(p.lambda_Y>=0 && p.zeta_Y0>=0, ...
    'Resistance-evolution parameters must be nonnegative.');
assert(p.R_H>=0 && p.R_tau0>=0 && p.R_dissipation>=0 && ...
    p.R_hardening>0 && ...
    p.R_max_fraction>0 && p.R_max_fraction<=1 && ...
    p.R_latent_scale>=0 && p.R_latent_scale<=1 && ...
    p.R_training_gain>=0 && p.R_training_rate>=0 && ...
    isfinite(p.R_entropy) && isfinite(p.R_Tref), ...
    'R-phase parameters are not admissible.');
assert(p.arp_tau0>=0 && p.arp_tau0_f>=0 && ...
    p.arp_tau0_f<=p.arp_tau0 && ...
    p.arp_tau_fraction_start>=0 && p.arp_tau_fraction_start<=1 && ...
    p.arp_tau_fraction_width>0 && p.arp_H>0 && ...
    p.arp_matrix_power>=0, ...
    'ARP parameters are not admissible.');
assert(p.mt_H>=0 && p.mt_hfor>=0 && ...
    p.mt_hfor_training_extra>=0 && p.mt_hfor_training_scale>0 && ...
    p.mt_hfor_training_power>=0 && ...
    p.mt_hfor_transition_dip>=0 && p.mt_hfor_transition_width>0 && ...
    p.mt_hfor_shoulder>=0 && p.mt_hfor_shoulder_start>=0 && ...
    p.mt_hfor_shoulder_start<=1 && p.mt_hfor_shoulder_width>0 && ...
    p.mt_hfor_shoulder_memory_scale>0 && ...
    p.mt_hfor_shoulder_memory_power>=0 && ...
    p.mt_hfor_shoulder_training_scale>0 && ...
    p.mt_hfor_shoulder_training_power>=0 && ...
    p.mt_hfor_activity>=0 && ...
    p.mt_hfor_activity_power>=1 && p.mt_hrev>=0 && ...
    p.mt_hrev_span>=0 && p.mt_hrev_span_load_power>=0 && ...
    p.mt_hrev_finish>=0 && p.mt_hrev_finish_low>=0 && ...
    p.mt_hrev_finish_low_load_power>=0 && ...
    p.mt_hrev_finish_low_power>=1 && ...
    p.mt_hrev_finish_low_training_scale>0 && ...
    p.mt_hrev_finish_low_training_power>=0 && ...
    p.mt_hrev_finish_load_power>=0 && ...
    p.mt_hrev_finish_power>=1 && ...
    p.mt_reverse_memory_shift>=0 && p.mt_reverse_memory_scale>0 && ...
    p.mt_reverse_memory_power>=0 && ...
    p.mt_n>=1, ...
    'MT nonlinear-hardening parameters are not admissible.');
assert(p.mt_Y_d>=0 && p.mt_Y_n>=1 && p.mt_Yp_sat>=0 && ...
    p.mt_Yp_high_sat>=0 && p.mt_Yp_memory_scale>0 && ...
    p.mt_Yp_memory_power>=0 && p.mt_Yp_activity_scale>0 && ...
    p.mt_Yp_activity_power>=0 && ...
    p.mt_Yp_rate>=0 && p.mt_Yp_reverse_scale>=0 && ...
    p.mt_Yp_reverse_training_decay>=0 && ...
    p.mt_Yp_reverse_memory_scale>0 && ...
    p.mt_Yp_reverse_memory_power>=0 && ...
    p.mt_training_reference_fraction>0 && ...
    p.mt_training_cycle_increment>=0 && ...
    p.mt_training_completion_tol>=p.numerics.xi_tol, ...
    'MT dissipation-evolution parameters are not admissible.');
assert(p.residual_xi_a>=0 && p.residual_xi_offset>=0 && ...
    p.residual_xi_peak_gain>=0 && ...
    p.residual_xi_offset+p.residual_xi_peak_gain<1, ...
    'Residual-martensite parameters are not admissible.');
end
