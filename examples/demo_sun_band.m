function geometry = demo_sun_band()
%DEMO_SUN_BAND Sample the elastic band shape factor and energy gradient.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'matlab'));
radius_m = 1e-4;
modulus_MPa = 60000;
poisson = 0.3;
radial_eigenstrain = 0.02;
length_ratio = [0.25; 0.5; 1; 2];
shape_factor = zeros(size(length_ratio));
gradient = zeros(size(length_ratio));
for k = 1:numel(length_ratio)
    shape_factor(k) = sma.candidate_sun_band_shape_factor_1d(length_ratio(k), poisson);
    gradient(k) = sma.candidate_sun_band_energy_gradient_1d( ...
        radius_m, modulus_MPa, poisson, radial_eigenstrain, length_ratio(k));
end
assert(all(isfinite(gradient)) && all(shape_factor > 0), 'The band calculation left its elastic domain.');
geometry = table(length_ratio, shape_factor, gradient, ...
    'VariableNames', {'LengthRadiusRatio', 'ShapeFactor', 'EnergyGradient_MPa_m2'});
disp(geometry);
end
