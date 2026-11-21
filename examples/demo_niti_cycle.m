function history = demo_niti_cycle(output_dir)
%DEMO_NITI_CYCLE Simulate a prescribed-strain NiTi cycle at fixed temperature.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'matlab'));
if nargin < 1
    output_dir = fullfile(root, 'outputs');
end

p = sma.params('niti_zhang2022');
state = sma.initial_state_1d(p);
temperature = state.T;
dt = 0.1;
loading = linspace(0, 0.04, 81).';
target = [loading; flipud(loading(1:end-1))];
time_s = (0:numel(target)-1).' * dt;
trajectory = zeros(numel(target), 4);
trajectory(1, :) = [state.sigma, state.eps, state.xi, state.T];
for k = 2:numel(target)
    trial = sma.update_1d_strain(p, state, target(k), temperature, dt, false);
    trajectory(k, :) = [trial.sigma, trial.eps, trial.xi, trial.T];
    state = trial;
end

assert(all(isfinite(trajectory(:))), 'The cycle contains a nonfinite state.');
assert(max(abs(trajectory(:, 2) - target)) <= 1e-8, 'Strain control did not close.');
assert(all(trajectory(:, 3) >= 0 & trajectory(:, 3) <= 1), 'Phase fraction is outside [0,1].');
assert(max(trajectory(:, 3)) > 0, 'The cycle did not activate martensite.');
history = table(time_s, trajectory(:, 1), trajectory(:, 2), trajectory(:, 3), trajectory(:, 4), ...
    'VariableNames', {'Time_s', 'Stress_MPa', 'LogarithmicStrain', 'MartensiteFraction', 'Temperature_K'});

if ~isfolder(output_dir)
    mkdir(output_dir);
end
writetable(history, fullfile(output_dir, 'niti-isothermal-cycle.csv'));
fig = figure('Visible', 'off', 'Color', 'white', 'Position', [100, 100, 760, 480]);
cleanup = onCleanup(@() close(fig));
ink = [0.12, 0.15, 0.19];
ax = axes(fig, 'Color', 'white', 'XColor', ink, 'YColor', ink, ...
    'GridColor', [0.74, 0.78, 0.82], 'GridAlpha', 0.6, ...
    'Box', 'off', 'TickDir', 'out', 'FontSize', 12);
hold(ax, 'on');
loading_rows = 1:numel(loading);
unloading_rows = numel(loading):height(history);
plot(ax, 100*history.LogarithmicStrain(loading_rows), history.Stress_MPa(loading_rows), ...
    'Color', [0.13, 0.34, 0.55], 'LineWidth', 2, 'DisplayName', 'Loading');
plot(ax, 100*history.LogarithmicStrain(unloading_rows), history.Stress_MPa(unloading_rows), ...
    'Color', [0.67, 0.31, 0.17], 'LineWidth', 2, 'DisplayName', 'Unloading');
xlabel(ax, 'Logarithmic strain (%)', 'Color', ink);
ylabel(ax, 'Stress (MPa)', 'Color', ink);
title(ax, sprintf('Isothermal NiTi cycle at %g K', temperature), 'Color', ink);
legend(ax, 'Location', 'northwest', 'TextColor', ink, 'Color', 'white', 'Box', 'off');
grid(ax, 'on');
exportgraphics(ax, fullfile(output_dir, 'niti-isothermal-cycle.png'), ...
    'Resolution', 130, 'BackgroundColor', 'white');
fprintf('NiTi | %d accepted increments | peak martensite fraction %.5f\n', numel(target)-1, max(history.MartensiteFraction));
end
