function factor = phase_activation_factor(p, trial_surface_overdrive_MPa)
%SMA.PHASE_ACTIVATION_FACTOR Report normalized phase-surface overdrive.

enabled = isfield(p, 'phase_overdrive_activation_enabled') && ...
    logical(p.phase_overdrive_activation_enabled);
if ~enabled
    factor = 1.0;
    return;
end

width_MPa = p.phase_overdrive_width_MPa;
assert(isfinite(width_MPa) && width_MPa > 0.0, ...
    'sma:phase_activation_factor:InvalidWidth', ...
    'The phase-overdrive activation width must be positive and finite.');
normalized = min(max(double(trial_surface_overdrive_MPa) / ...
    width_MPa, 0.0), 1.0);
shape = 'cubic_smoothstep';
if isfield(p, 'phase_overdrive_activation_shape') && ...
        ~isempty(p.phase_overdrive_activation_shape)
    shape = char(string(p.phase_overdrive_activation_shape));
end
switch lower(shape)
    case 'cubic_smoothstep'
        factor = normalized.^2 .* (3.0 - 2.0 .* normalized);
    case 'linear_minimax_ramp'
        factor = normalized;
    otherwise
        error('sma:phase_activation_factor:UnknownShape', ...
            'Unsupported phase-overdrive activation shape: %s.', shape);
end
end
