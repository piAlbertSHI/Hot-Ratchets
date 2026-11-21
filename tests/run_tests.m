function run_tests()
%RUN_TESTS Check the public material kernels and synthetic histories.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'matlab'), fullfile(root, 'tests'));
checks = {@verify_force_history_step_contract, @verify_strain_history_step_contract, ...
    @verify_mixed_mt_properties_contract, @verify_unified_state_contract_1d, ...
    @verify_restart_schema_v1, @verify_gradient_stress_thermomechanical_control_contract, ...
    @verify_yuan_reduced_material_provider_contract, @verify_sun_band_shape_factor_contract, ...
    @verify_sun_band_energy_gradient_contract};
for k = 1:numel(checks)
    feval(checks{k});
    fprintf('PASS %s\n', func2str(checks{k}));
end
fprintf('%d public contract checks passed.\n', numel(checks));
end
