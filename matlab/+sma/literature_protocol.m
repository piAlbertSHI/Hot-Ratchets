function protocol = literature_protocol(case_name, run_mode, parameter_overrides)
%SMA.LITERATURE_PROTOCOL Published loading histories for short or long runs.
% Modes: 0 stress, 1 strain, 2 stress-rate unload, 4 strain envelope.

if nargin < 2
    run_mode = 'short';
end
if nargin < 3
    parameter_overrides = struct();
end
case_id = sma.case_id(case_name);
p = sma.params_from_id(case_id);
p = apply_overrides(p,parameter_overrides);
is_long = strcmpi(run_mode, 'long');

switch case_id
    case 1
        % The published fast-unloading path contains several intermediate
        % points. A 0.25 s top-level step preserves the 80 MPa/s rate while
        % resolving the plotted reverse-transformation branch. Zhang reports
        % total unload times of 4, 6, 8, 10 and 12 s for S4--S12; the lower
        % bound is therefore 4 s rather than the former unsupported 8 s.
        dt = 0.25;
        n_cycles = choose(is_long, 20, 3);
        peak_strain = niti_peak_strain(case_name);
        load_time_s = 1.0e3*peak_strain;
        unload_time_s = max(1.0e2*peak_strain,4.0);
        base = [4 0.000 peak_strain load_time_s; ...
            1 peak_strain peak_strain 60.0; ...
            2 80.000 80.000 unload_time_s; ...
            0 0.000 0.000 60.0];
        label = sprintf('Zhang 2022 Figure 2 S%d', ...
            round(100*peak_strain));
    case {2,3,4}
        dt = 0.25;
        if case_id == 2
            peak = 150.0;
            n_cycles = choose(is_long, 50, 3);
        elseif case_id == 3
            peak = 300.0;
            n_cycles = choose(is_long, 100, 3);
        else
            peak = 300.0;
            n_cycles = choose(is_long, 275, 3);
        end
        base = [0 0.000 peak 1.0; 0 peak peak 15.0; ...
            0 peak 0.000 1.0; 0 0.000 0.000 15.0];
        label = 'Yuan 2020 Figure 6';
    case 5
        dt = 0.50;
        n_cycles = choose(is_long, 100, 3);
        base = [1 0.000 0.040 1.5; 1 0.040 0.040 28.5; ...
            2 80.000 80.000 1.5; 0 0.000 0.000 28.5];
        label = ['Imran and Zhang 2020 Figure 5(d,e), ' ...
            'secondary-source trend benchmark'];
    case 7
        dt = 0.50;
        n_cycles = choose(is_long, 200, 3);
        base = [0 0.000 65.000 1.5; 0 65.000 65.000 28.5; ...
            2 (65/1.5) (65/1.5) 1.5; 0 0.000 0.000 28.5];
        label = ['Imran and Zhang 2020 Figure 4, ' ...
            'secondary-source trend benchmark'];
    otherwise
        error('sma:literature_protocol:Unsupported', ...
            'No literature protocol for case ID %d.', case_id);
end

segments = make_segments(base, n_cycles);
h = compile_segments(segments, dt);
protocol = h;
protocol.case_id = case_id;
protocol.case_name = char(case_name);
protocol.label = label;
protocol.p = p;
protocol.dt = dt;
protocol.n_cycles = n_cycles;
protocol.run_mode = char(run_mode);
protocol.stop_time = h.t(end);
protocol.signal = timeseries([h.mode,h.target],h.t);
% Keep source facts and validation eligibility attached to every protocol
% object.  Numerical runners may still execute the historical surrogate, but
% downstream reports can no longer mistake it for a source-authoritative
% experimental history.
protocol.authority = sma.protocol_authority(case_name);
% Keep source and execution strain coordinates explicit in the protocol
% transaction.  The historical numeric target remains untouched; the
% source-declared branch is an opt-in mapping for Zhang's current-length
% strain x=(l-l0)/l to the additive engineering target e=x/(1-x).
protocol.source_target = protocol.target;
protocol.execution_target_legacy = protocol.target;
protocol.execution_target_source_declared = protocol.target;
protocol.execution_strain_measure_legacy = 'additive_small_strain_work_conjugate';
protocol.execution_strain_measure_source_declared = 'engineering_reference_strain';
protocol.source_to_execution_strain_map = 'identity_numeric_default';
if case_id==1 && strcmp(protocol.authority.source_strain_measure, ...
        'current_length_normalized_(l-l0)/l')
    mask=ismember(protocol.mode,[1,4]);
    x=protocol.target(mask);
    assert(all(x>-1 & x<1), ...
        'Zhang source current-length target must remain in (-1,1).');
    protocol.execution_target_source_declared(mask)=x./(1-x);
    protocol.source_to_execution_strain_map= ...
        'current_length_to_engineering_e_equals_x_over_1_minus_x';
end
end

function p = apply_overrides(p,overrides)
names = fieldnames(overrides);
for k = 1:numel(names)
    name = names{k};
    assert(isfield(p,name), ...
        'sma:literature_protocol:UnknownParameter', ...
        'Unknown material parameter override: %s.',name);
    p.(name) = overrides.(name);
end
assert(p.EA>0 && p.EM>0 && p.cp>0 && p.rho>0, ...
    'Material moduli, density and heat capacity must be positive.');
assert(p.Hmax_i>=p.Hmin && p.Hmax_f>=p.Hmin, ...
    'Hmax must not be below Hmin.');
assert(p.latent>0 && p.latent_f>0 && p.latent_f<=p.latent, ...
    'Latent-heat overrides are not admissible.');
end

function value = niti_peak_strain(case_name)
name = lower(char(case_name));
tokens = regexp(name,'niti_s(4|6|8|10|12)$','tokens','once');
if isempty(tokens)
    value = 0.08;
else
    value = str2double(tokens{1})/100.0;
end
end

function value = choose(condition, true_value, false_value)
if condition
    value = true_value;
else
    value = false_value;
end
end

function segments = make_segments(base, n_cycles)
n_stages = size(base,1);
segments = repmat([base,zeros(n_stages,2)],n_cycles,1);
segments(:,5) = repelem((1:n_cycles).',n_stages);
segments(:,6) = repmat((1:n_stages).',n_cycles,1);
end

function h = compile_segments(segments, dt)
steps = max(1,round(segments(:,4)/dt));
n_total = sum(steps);
h = struct('t',zeros(n_total+1,1),'mode',zeros(n_total+1,1), ...
    'target',zeros(n_total+1,1),'cycle',zeros(n_total+1,1), ...
    'stage',zeros(n_total+1,1));
h.cycle(1) = 1;
cursor = 1;
for k = 1:size(segments,1)
    n = steps(k);
    indices = cursor + (1:n);
    values = linspace(segments(k,2),segments(k,3),n+1).';
    h.t(indices) = h.t(cursor) + (1:n).' * dt;
    h.mode(indices) = segments(k,1);
    h.target(indices) = values(2:end);
    h.cycle(indices) = segments(k,5);
    h.stage(indices) = segments(k,6);
    cursor = cursor + n;
end
end
