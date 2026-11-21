function [trial_matrix,force_N,report] = update_series_nodes( ...
    p,committed_matrix,temperature_K,length_m,area_m2, ...
    target_extension_m,dt,force_guess_N)
%SMA.UPDATE_SERIES_NODES Rollback-safe common-force compatibility Newton.
%
% Every trial force is evaluated from COMMITTED_MATRIX. No node state is
% advanced until the common-force residual converges.

node_count = size(committed_matrix,1);
state_width=sma.state_width_1d(p);
assert(size(committed_matrix,2)==state_width, ...
    'State matrix must be N-by-%d.',state_width);
temperature_K = expand_vector(temperature_K,node_count,'temperature');
length_m = expand_vector(length_m,node_count,'length');
area_m2 = expand_vector(area_m2,node_count,'area');
assert(all(area_m2>0) && all(length_m>0),'Geometry must be positive.');
if nargin < 8 || isempty(force_guess_N)
    stress_mean_MPa = mean(committed_matrix(:,7));
    force_guess_N = stress_mean_MPa*1e6/mean(1.0./area_m2);
end

force_N = force_guess_N;
converged = false;
history = zeros(p.numerics.max_iterations,3);
for iteration = 1:p.numerics.max_iterations
    [trial_matrix,residual_m,derivative_m_per_N,max_local_residual] = ...
        evaluate_force(p,committed_matrix,temperature_K,length_m, ...
        area_m2,target_extension_m,dt,force_N);
    history(iteration,:) = [force_N,residual_m,max_local_residual];
    if abs(residual_m)<=1e-10 && max_local_residual<5e-6
        converged = true;
        break;
    end
    if ~isfinite(derivative_m_per_N) || derivative_m_per_N<=0
        error('sma:update_series_nodes:InvalidTangent', ...
            'Compatibility tangent must be positive.');
    end
    step_N = -residual_m/derivative_m_per_N;
    eta = 1.0;
    accepted = false;
    old_norm = abs(residual_m);
    for ls = 1:p.numerics.max_line_search
        candidate_force = force_N+eta*step_N;
        [~,candidate_residual] = evaluate_force( ...
            p,committed_matrix,temperature_K,length_m,area_m2, ...
            target_extension_m,dt,candidate_force);
        if abs(candidate_residual)<old_norm
            force_N = candidate_force;
            accepted = true;
            break;
        end
        eta = 0.5*eta;
    end
    if ~accepted
        break;
    end
end
if ~converged
    error('sma:update_series_nodes:Compatibility', ...
        'Series compatibility failed with residual %.6g m.',residual_m);
end
history = history(1:iteration,:);
report = struct( ...
    'converged',converged, ...
    'iterations',iteration, ...
    'compatibility_residual_m',residual_m, ...
    'max_local_residual',max_local_residual, ...
    'force_history',history, ...
    'committed_state_unchanged',true);
end

function [trial_matrix,residual_m,derivative_m_per_N,max_local_residual] = ...
    evaluate_force(p,committed_matrix,temperature_K,length_m,area_m2, ...
    target_extension_m,dt,force_N)
node_count = size(committed_matrix,1);
trial_matrix = zeros(size(committed_matrix));
extension_m = 0.0;
derivative_m_per_N = 0.0;
max_local_residual = 0.0;
for i = 1:node_count
    committed = sma.unpack_state_1d(p,committed_matrix(i,:).');
    sigma_MPa = force_N/area_m2(i)/1e6;
    [trial,aux] = sma.update_1d_stress( ...
        p,committed,sigma_MPa,temperature_K(i),dt,true);
    trial_matrix(i,:) = sma.pack_state_1d(p,trial).';
    extension_m = extension_m+length_m(i)*trial.eps;
    derivative_m_per_N = derivative_m_per_N+ ...
        length_m(i)*aux.deps_dsigma/(area_m2(i)*1e6);
    max_local_residual = max(max_local_residual,aux.residual_norm);
end
residual_m = extension_m-target_extension_m;
end

function value = expand_vector(value,node_count,label)
value = value(:);
if isscalar(value)
    value = repmat(value,node_count,1);
end
assert(numel(value)==node_count,'%s width mismatch.',label);
end
