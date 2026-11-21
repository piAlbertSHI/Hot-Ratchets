function history = demo_cualmn_ramp()
%DEMO_CUALMN_RAMP Advance the reduced Cu-Al-Mn material under a strain ramp.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'matlab'));
p = sma.params('cualmn_annealed');
state = sma.initial_state_1d(p);
dt = 0.01;
control = struct('mode', 'engineering_strain_rate', 'signed_rate_1_per_s', 0.15);
increments = 4;
values = zeros(increments + 1, 4);
values(1, :) = [expm1(state.eps), state.sigma, state.xi, state.T];
for k = 1:increments
    [trial, receipt] = sma.candidate_yuan_reduced_material_step_1d(p, state, control, p.Tamb, dt);
    assert(receipt.accepted, 'The reduced Cu-Al-Mn step did not close.');
    values(k+1, :) = [expm1(trial.eps), trial.sigma, trial.xi, trial.T];
    state = trial;
end
assert(all(isfinite(values(:))), 'The ramp contains a nonfinite state.');
assert(abs(values(end, 1) - increments*dt*control.signed_rate_1_per_s) <= 1e-9, ...
    'Engineering strain control did not close.');
history = table((0:increments).'*dt, values(:, 1), values(:, 2), values(:, 3), values(:, 4), ...
    'VariableNames', {'Time_s', 'EngineeringStrain', 'Stress_MPa', 'ReactionExtent', 'Temperature_K'});
disp(history);
end
