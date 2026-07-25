function [trial, aux] = update_3d_strain(p, state, h_target, T, dt)
%SMA.UPDATE_3D_STRAIN Strain-driven wrapper around the unique tensor core.
% The outer stress solve always evaluates trials from the same committed state.

if nargin < 5 || isempty(dt)
    dt = p.default_dt;
end
target = sma.tensor_to_voigt(h_target);
tau0 = sma.tensor_to_voigt(state.tau);
fun = @(tau_v) strain_residual(p, state, tau_v, target, T, dt);
tau_v = tau0;
converged = false;
for iteration = 1:p.numerics.max_iterations
    R = fun(tau_v);
    if norm(R, inf) <= p.numerics.residual_tol
        converged = true;
        break;
    end
    J = sma.numerical_jacobian(fun, tau_v, p.numerics.jacobian_step);
    step = -J \ R;
    eta = 1.0;
    old_norm = norm(R, inf);
    accepted = false;
    for ls = 1:p.numerics.max_line_search
        candidate = tau_v + eta * step;
        if norm(fun(candidate), inf) < old_norm
            tau_v = candidate;
            accepted = true;
            break;
        end
        eta = 0.5 * eta;
    end
    if ~accepted
        break;
    end
end
if ~converged
    error('sma:update_3d_strain:OuterNewton', ...
        'Strain-driven stress solve failed with residual %g.', norm(fun(tau_v), inf));
end
[trial, aux] = sma.update_3d_stress( ...
    p, state, sma.voigt_to_tensor(tau_v), T, dt, true);
aux.outer_iterations = iteration;
aux.strain_residual = sma.tensor_to_voigt(trial.h) - target;
if all(isfinite(aux.compliance_alg), 'all')
    aux.tangent_alg = inv(aux.compliance_alg);
else
    aux.tangent_alg = nan(6);
end
end

function R = strain_residual(p, state, tau_v, target, T, dt)
[trial, ~] = sma.update_3d_stress( ...
    p, state, sma.voigt_to_tensor(tau_v), T, dt, false);
R = sma.tensor_to_voigt(trial.h) - target;
end
