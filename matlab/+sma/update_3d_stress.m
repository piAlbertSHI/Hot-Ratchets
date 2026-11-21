function [trial, aux] = update_3d_stress(p, state, tau, T, dt, compute_tangent)
%SMA.UPDATE_3D_STRESS Tensor mother-model update at prescribed Kirchhoff stress.

if nargin < 6
    compute_tangent = false;
end
if nargin < 5 || isempty(dt)
    dt = p.default_dt;
end
dt = max(dt, eps);
assert(isequal(size(tau), [3, 3]) && norm(tau - tau.', 'fro') < 1e-10, ...
    'Kirchhoff stress must be a symmetric 3x3 tensor.');
state=sma.ensure_population_state(p,state);
state=sma.ensure_zhang_direct_state(p,state);
state=sma.ensure_r_phase_state(p,state);
state=canonicalize_derived_state(p,state);

[seq_trial, ~] = sma.equivalent_stress(tau + state.beta);
[phi_f, phi_r] = sma.trial_surfaces_3d(p,state,tau,T);
branch = choose_branch(p, state, seq_trial, T, phi_f, phi_r);

if branch == 0
    trial = elastic_state(p, state, tau, T);
    aux = make_aux(p, state, trial, dt, true, 0, 0.0, 0.0);
    aux.kinetic_rate_1_per_s = 0.0;
    aux.kinetic_fraction = 0.0;
    aux.thermodynamic_target_xi = state.xi;
    aux.p_AM_coefficient_1_per_s = 0.0;
    aux.p_MA_coefficient_1_per_s = 0.0;
    aux.dxi_dt_1_per_s = 0.0;
    aux.kinetic_branch = 0.0;
    aux.phase_activation_factor = 0.0;
    aux.trial_surface_overdrive_MPa = max([phi_f, phi_r, 0.0]);
else
    reverse = reverse_reference(state, branch);
    [x0, ~] = initial_guess(p, state, tau, T, branch);
    if p.enable_composition_weighted_mt
        x0(1)=active_surface_fraction_seed( ...
            p,state,tau,T,branch,reverse,x0(1));
    end
    boundary_probe = x0;
    boundary_probe(1) = double(branch > 0);
    boundary_R = sma.residual_3d( ...
        p,state,tau,T,branch,reverse,boundary_probe);
    boundary = boundary_R(1) > p.numerics.surface_tol;
    if boundary
        [x, converged, iterations, residual_norm] = solve_boundary( ...
            p, state, tau, T, branch, reverse, x0);
    else
        [x, converged, iterations, residual_norm] = solve_active( ...
            p, state, tau, T, branch, reverse, x0);
        if ~converged
            [boundary_x, boundary_converged, boundary_iterations, ...
                boundary_residual] = ...
                solve_boundary(p, state, tau, T, branch, reverse, x);
            boundary_R = sma.residual_3d( ...
                p,state,tau,T,branch,reverse,boundary_x);
            boundary_admissible = boundary_R(1) >= ...
                -p.numerics.surface_tol;
            if boundary_converged && boundary_admissible
                x = boundary_x;
                converged = true;
                residual_norm = boundary_residual;
                iterations = iterations + boundary_iterations;
            end
        end
    end
    if ~converged && p.numerics.cutback_on_failure
        error('sma:update_3d_stress:LocalNewton', ...
            '3D local Newton failed: branch=%d, residual=%g.', branch, residual_norm);
    end
    [~, data] = sma.residual_3d(p, state, tau, T, branch, reverse, x);
    trial = state;
    trial.T = T;
    trial.xi = min(max(x(1), 0.0), 1.0);
    trial.h_tr = sma.voigt_to_tensor(x(2:7));
    trial.h_tp = sma.voigt_to_tensor(x(8:13));
    trial.zeta_d = max(x(14), state.zeta_d);
    trial.beta = sma.voigt_to_tensor(x(15:20));
    trial.tau = tau;
    trial.branch = branch;
    trial.N = data.N;
    trial.xi_d = data.xi_d;
    trial.Hcur = data.Hcur;
    trial.Hmax = data.Hmax;
    trial.xi_res = state.xi_res;
    if branch > 0
        trial.xi_r = trial.xi;
        trial.h_tr_r = trial.h_tr;
    else
        trial.xi_res = sma.residual_martensite_limit(p,state);
        if trial.xi <= trial.xi_res+p.numerics.xi_tol
            trial.xi = trial.xi_res;
        end
    end
    trial.Y_population_activity=sma.update_population_activity( ...
        p,state.Y_population_activity,state.xi,trial.xi);
    history_final=sma.candidate_zhang_history_3d( ...
        p,state,trial.xi,trial.h_tp,tau);
    trial=sma.commit_zhang_direct_history(trial,history_final);
    trial.h_arp = sma.update_arp_3d(p,state,tau,trial.xi);
    trial.zeta_p = sma.update_training_activity(p,state,trial,branch);
    thermodynamic_target_xi = trial.xi;
    trial_surface_overdrive_MPa = choose(branch > 0, phi_f, phi_r);
    [trial, kinetic_rate, kinetic_fraction, activation_factor] = ...
        apply_rate_regularization_3d( ...
        p,state,trial,branch,tau,dt,trial_surface_overdrive_MPa);
    alpha = sma.thermal_expansion(p,trial.xi);
    [h_R,xi_R_trial] = sma.r_phase_strain_3d( ...
        p,tau,T,trial.zeta_p,trial.xi,state.xi_R,state.xi);
    if p.enable_explicit_r_phase_state, trial.xi_R=xi_R_trial; end
    trial.h = sma.elastic_strain_3d(p,trial.xi,tau)+ ...
        alpha*(T-p.T0)*eye(3)+h_R+trial.h_tr+trial.h_tp+trial.h_arp;
    aux = make_aux(p,state,trial,dt,converged,iterations, ...
        residual_norm,data.phi);
    aux.kinetic_rate_1_per_s = kinetic_rate;
    aux.kinetic_fraction = kinetic_fraction;
    aux.thermodynamic_target_xi = thermodynamic_target_xi;
    if isfinite(kinetic_rate)
        aux.p_AM_coefficient_1_per_s = ...
            kinetic_rate*thermodynamic_target_xi;
        aux.p_MA_coefficient_1_per_s = ...
            kinetic_rate*(1.0-thermodynamic_target_xi);
    else
        aux.p_AM_coefficient_1_per_s = 0.0;
        aux.p_MA_coefficient_1_per_s = 0.0;
    end
    aux.dxi_dt_1_per_s = (trial.xi-state.xi)/dt;
    aux.kinetic_branch = branch;
    aux.phase_activation_factor = activation_factor;
    aux.trial_surface_overdrive_MPa = trial_surface_overdrive_MPa;
end

if compute_tangent
    aux.compliance_alg = numerical_compliance(p, state, tau, T, dt);
else
    aux.compliance_alg = nan(6);
end

function [trial, rate_1_per_s, kinetic_fraction, activation_factor] = ...
        apply_rate_regularization_3d( ...
        p,state,trial,branch,tau,dt,trial_surface_overdrive_MPa)
rate_1_per_s = inf;
kinetic_fraction = 1.0;
activation_factor = 1.0;
if ~p.enable_rate_regularization
    return;
end
if branch > 0
    base_rate = p.phase_loading_rate_1_per_s;
else
    base_rate = p.phase_unloading_rate_1_per_s;
end
assert(isfinite(base_rate) && base_rate > 0.0, ...
    'sma:update_3d_stress:InvalidRateRegularization', ...
    'The selected finite directional rate must be positive and finite.');
activation_factor = sma.phase_activation_factor( ...
    p, trial_surface_overdrive_MPa);
assert(isfield(p, 'phase_rate_multiplier_mode') && strcmp( ...
    char(string(p.phase_rate_multiplier_mode)), ...
    'none_qian_eq3_direction_coefficient_exact'), ...
    'sma:update_3d_stress:UnsupportedPhaseRateMultiplier', ...
    'The unique Qian Eq. (3) model requires an exact directional rate.');
% Keep the tensor mother model identical to the 1D projection: branch and
% target come from the return map, while Qian Eq. (3) retains its fixed
% directional coefficient independently of surface distance.
rate_1_per_s = base_rate;
kinetic_fraction = -expm1(-rate_1_per_s*dt);
kinetic_fraction = min(max(kinetic_fraction, 0.0), 1.0);
if kinetic_fraction >= 1.0-10.0*eps
    return;
end

target = trial;
trial.xi = state.xi + kinetic_fraction*(target.xi-state.xi);
trial.h_tr = state.h_tr + kinetic_fraction*(target.h_tr-state.h_tr);
trial.h_tp = state.h_tp + kinetic_fraction*(target.h_tp-state.h_tp);
trial.zeta_d = max(state.zeta_d, state.zeta_d + ...
    kinetic_fraction*(target.zeta_d-state.zeta_d));
trial.beta = state.beta + kinetic_fraction*(target.beta-state.beta);
trial.xi_res = target.xi_res;
if branch > 0
    trial.xi_r = trial.xi;
    trial.h_tr_r = trial.h_tr;
else
    trial.xi_r = state.xi_r;
    trial.h_tr_r = state.h_tr_r;
end
trial.Y_population_activity=sma.update_population_activity( ...
    p,state.Y_population_activity,state.xi,trial.xi);
history_final=sma.candidate_zhang_history_3d( ...
    p,state,trial.xi,trial.h_tp,tau);
trial=sma.commit_zhang_direct_history(trial,history_final);
trial.h_arp = sma.update_arp_3d(p,state,tau,trial.xi);
trial.zeta_p = sma.update_training_activity(p,state,trial,branch);
[~,xi_R_trial] = sma.r_phase_strain_3d( ...
    p,tau,T,trial.zeta_p,trial.xi,state.xi_R,state.xi);
if p.enable_explicit_r_phase_state, trial.xi_R=xi_R_trial; end
[effective_stress_MPa, trial.N] = ...
    sma.equivalent_stress(tau+trial.beta);
[trial.Hcur,trial.Hmax] = sma.hcur( ...
    p,effective_stress_MPa,trial.zeta_d);
trial.xi_d = min(max(trial.Hcur/max(trial.Hmax,eps),0.0),1.0)*trial.xi;
end
end

function branch = choose_branch(p, state, seq, T, phi_f, phi_r)
tol = p.numerics.surface_tol;
[old_seq, ~] = sma.equivalent_stress(state.tau + state.beta);
loading_measure = seq-old_seq;
temperature_change = T-state.T;
can_forward = state.xi < 1.0 - p.numerics.xi_tol && phi_f > tol;
reverse_limit = sma.residual_martensite_limit(p,state);
can_reverse = state.xi > reverse_limit+p.numerics.xi_tol && phi_r > tol;
if p.enable_surface_driven_branch_selection
    if can_forward && can_reverse
        branch=choose(phi_f>=phi_r,1,-1);
    elseif can_forward
        branch=1;
    elseif can_reverse
        branch=-1;
    else
        branch=0;
    end
    return;
end
mechanical_tol = p.numerics.direction_tol;
if loading_measure > mechanical_tol
    if can_forward
        branch = 1;
    else
        branch = 0;
    end
elseif loading_measure < -mechanical_tol
    if can_reverse
        branch = -1;
    else
        branch = 0;
    end
else
    if can_forward && can_reverse
        branch = choose(temperature_change<0,1,-1);
    elseif can_forward
        branch = 1;
    elseif can_reverse
        branch = -1;
    else
        branch = 0;
    end
end
end

function value = choose(condition,true_value,false_value)
if condition
    value = true_value;
else
    value = false_value;
end
end

function reverse = reverse_reference(state, branch)
reverse = struct('xi_r', state.xi_r, 'h_tr_r', state.h_tr_r);
if branch < 0 && reverse.xi_r <= 0.0
    reverse.xi_r = state.xi;
    reverse.h_tr_r = state.h_tr;
end
end

function [x0, boundary] = initial_guess(p, state, tau, T, branch)
[seq, ~] = sma.equivalent_stress(tau + state.beta);
if branch > 0
    xi_target = (seq / p.CM - (T - p.Ms)) / (p.Ms - p.Mf);
    xi_target = max(xi_target, state.xi);
    boundary = xi_target >= 1.0;
else
    xi_target = (seq / p.CA - (T - p.Af)) / (p.Af - p.As);
    xi_target = min(xi_target, state.xi);
    reverse_limit = sma.residual_martensite_limit(p,state);
    boundary = xi_target <= reverse_limit;
    xi_target = max(xi_target,reverse_limit);
end
xi_target = min(max(xi_target, 0.0), 1.0);
x0 = [xi_target; sma.tensor_to_voigt(state.h_tr); ...
    sma.tensor_to_voigt(state.h_tp); state.zeta_d; ...
    sma.tensor_to_voigt(state.beta)];
end

function xi_seed=active_surface_fraction_seed( ...
    p,state,tau,T,branch,reverse,xi_fallback)
if branch>0
    lo=state.xi; hi=1.0;
else
    lo=sma.residual_martensite_limit(p,state); hi=state.xi;
end
surface=@(xi)surface_at_fraction(p,state,tau,T,branch,reverse,xi);
flo=surface(lo); fhi=surface(hi);
if ~isfinite(flo) || ~isfinite(fhi) || flo*fhi>0.0
    xi_seed=xi_fallback;
    return;
end
for k=1:48
    mid=0.5*(lo+hi); fmid=surface(mid);
    if abs(fmid)<=0.1*p.numerics.residual_tol
        lo=mid; hi=mid;
        break;
    elseif flo*fmid<=0.0
        hi=mid; fhi=fmid; %#ok<NASGU>
    else
        lo=mid; flo=fmid;
    end
end
xi_seed=0.5*(lo+hi);
end

function value=surface_at_fraction(p,state,tau,T,branch,reverse,xi)
x=[xi;sma.tensor_to_voigt(state.h_tr); ...
    sma.tensor_to_voigt(state.h_tp);state.zeta_d; ...
    sma.tensor_to_voigt(state.beta)];
R=sma.residual_3d(p,state,tau,T,branch,reverse,x);
value=R(1);
end

function [x, converged, iterations, residual_norm] = solve_active( ...
    p, state, tau, T, branch, reverse, x)
y = [x(1);state.zeta_d;sma.tensor_to_voigt(state.beta)];
fun = @(memory) active_memory_residual( ...
    p,state,tau,T,branch,reverse,memory);
converged = false;
for iterations = 1:p.numerics.max_iterations
    R = fun(y);
    residual_norm = norm(R, inf);
    if residual_norm <= p.numerics.residual_tol
        converged = true;
        break;
    end
    if branch > 0
        xi_lower = state.xi;
        xi_upper = 1.0;
    else
        xi_lower = sma.residual_martensite_limit(p,state);
        xi_upper = state.xi;
    end
    lower = [xi_lower;state.zeta_d;-inf(6,1)];
    upper = [xi_upper;inf;inf(6,1)];
    J = sma.box_constrained_jacobian( ...
        fun,y,p.numerics.jacobian_step,lower,upper);
    step = -J \ R;
    [y, accepted] = line_search_active( ...
        fun,y,step,state,branch,residual_norm,p);
    if ~accepted || norm(step, inf) <= p.numerics.step_tol
        break;
    end
end
iterations = min(iterations, p.numerics.max_iterations);
x_probe = [y(1);sma.tensor_to_voigt(state.h_tr); ...
    sma.tensor_to_voigt(state.h_tp);y(2);y(3:8)];
[~,data] = sma.residual_3d(p,state,tau,T,branch,reverse,x_probe);
delta_xi = y(1)-state.xi;
x = [y(1);sma.tensor_to_voigt(state.h_tr+data.Lambda*delta_xi); ...
    sma.tensor_to_voigt(state.h_tp+data.Lambda_tp*delta_xi); ...
    y(2);y(3:8)];
final_R = sma.residual_3d(p,state,tau,T,branch,reverse,x);
residual_norm = norm(final_R,inf);
converged = converged || residual_norm<=p.numerics.residual_tol;
end

function [x, converged, iterations, residual_norm] = solve_boundary( ...
    p, state, tau, T, branch, reverse, ~)
if branch>0
    fixed_xi = 1.0;
else
    fixed_xi = sma.residual_martensite_limit(p,state);
end
y = [state.zeta_d;sma.tensor_to_voigt(state.beta)];
converged = false;
for iterations = 1:100
    x_probe = [fixed_xi;sma.tensor_to_voigt(state.h_tr); ...
        sma.tensor_to_voigt(state.h_tp);y(1);y(2:7)];
    [~,data] = sma.residual_3d(p,state,tau,T,branch,reverse,x_probe);
    target = [state.zeta_d+abs(data.xi_d-state.xi_d); ...
        sma.tensor_to_voigt(data.beta_target)];
    residual_norm = norm(y-target,inf);
    if residual_norm<=p.numerics.residual_tol
        converged = true;
        break;
    end
    y = 0.5*y+0.5*target;
    y(1) = max(y(1),state.zeta_d);
end
x_probe = [fixed_xi;sma.tensor_to_voigt(state.h_tr); ...
    sma.tensor_to_voigt(state.h_tp);y(1);y(2:7)];
[~,data] = sma.residual_3d(p,state,tau,T,branch,reverse,x_probe);
delta_xi = fixed_xi-state.xi;
x = [fixed_xi; ...
    sma.tensor_to_voigt(state.h_tr+data.Lambda*delta_xi); ...
    sma.tensor_to_voigt(state.h_tp+data.Lambda_tp*delta_xi); ...
    y(1);y(2:7)];
final_R = sma.residual_3d(p,state,tau,T,branch,reverse,x);
if branch<0
    residual_norm = norm(final_R(8:20),inf);
else
    residual_norm = norm(final_R(2:20),inf);
end
converged = converged || residual_norm<=p.numerics.residual_tol;
end

function R = active_memory_residual( ...
    p,state,tau,T,branch,reverse,memory)
x = [memory(1);sma.tensor_to_voigt(state.h_tr); ...
    sma.tensor_to_voigt(state.h_tp);memory(2);memory(3:8)];
all_R = sma.residual_3d(p,state,tau,T,branch,reverse,x);
R = all_R([1,14:20]);
end

function [candidate, accepted] = line_search_active(fun, x, step, state, branch, ...
    residual_norm, p)
eta = 1.0;
accepted = false;
candidate = x;
for ls = 1:p.numerics.max_line_search
    value = x + eta * step;
    value(1) = min(max(value(1), 0.0), 1.0);
    value(2) = max(value(2), state.zeta_d);
    if branch > 0
        value(1) = max(value(1), state.xi);
    else
        reverse_limit = sma.residual_martensite_limit(p,state);
        value(1) = min(max(value(1),reverse_limit),state.xi);
    end
    if norm(fun(value), inf) < residual_norm
        candidate = value;
        accepted = true;
        return;
    end
    eta = 0.5 * eta;
end
end

function trial = elastic_state(p, state, tau, T)
trial = state;
trial.T = T;
trial.tau = tau;
trial.branch = 0.0;
[seq,~]=sma.equivalent_stress(tau+state.beta);
[trial.Hcur,trial.Hmax]=sma.hcur(p,seq,state.zeta_d);
trial.xi_d=min(max(trial.Hcur/max(trial.Hmax,eps),0.0),1.0)*state.xi;
trial.h_arp = sma.update_arp_3d(p,state,tau,state.xi);
alpha = sma.thermal_expansion(p, state.xi);
[h_R,xi_R_trial] = sma.r_phase_strain_3d( ...
    p,tau,T,state.zeta_p,state.xi,state.xi_R,state.xi);
if p.enable_explicit_r_phase_state, trial.xi_R=xi_R_trial; end
trial.h = sma.elastic_strain_3d(p, state.xi, tau) + ...
    alpha * (T - p.T0) * eye(3) + h_R + state.h_tr + ...
    state.h_tp + trial.h_arp;
end

function state=canonicalize_derived_state(p,state)
[seq,~]=sma.equivalent_stress(state.tau+state.beta);
[state.Hcur,state.Hmax]=sma.hcur(p,seq,state.zeta_d);
state.xi_d=min(max(state.Hcur/max(state.Hmax,eps),0.0),1.0)*state.xi;
end

function aux = make_aux(p, state, trial, dt, converged, iterations, ...
    residual_norm, phi)
delta_xi = trial.xi - state.xi;
delta_h_tr = trial.h_tr - state.h_tr;
delta_h_tp = trial.h_tp - state.h_tp;
delta_h_arp = trial.h_arp - state.h_arp;
latent_value = sma.effective_latent_heat(p,state,trial);
latent_specific = latent_value * delta_xi;
[seq_old,~] = sma.equivalent_stress(state.tau);
[seq_new,~] = sma.equivalent_stress(trial.tau);
r_latent_specific = sma.r_phase_latent_specific( ...
    p,state.xi,seq_old,state.T,state.zeta_p, ...
    trial.xi,seq_new,trial.T,trial.zeta_p,state.xi_R,trial.xi_R);
latent_specific = latent_specific+r_latent_specific;
work_mpa = sum((trial.tau + trial.beta) .* delta_h_tr, 'all') + ...
    sum(trial.tau .* (delta_h_tp+delta_h_arp), 'all');
mechanical_specific = 1.0e6 / p.rho * abs(work_mpa);
dissipation_specific = p.taylor_quinney * mechanical_specific;
strain_work_MPa = 0.5*sum( ...
    (state.tau+trial.tau).*(trial.h-state.h),'all');
aux = struct( ...
    'converged', converged, ...
    'iterations', iterations, ...
    'residual_norm', residual_norm, ...
    'phi', phi, ...
    'branch', trial.branch, ...
    'delta_xi', delta_xi, ...
    'effective_latent_Jkg', latent_value, ...
    'latent_specific_Jkg', latent_specific, ...
    'r_phase_latent_specific_Jkg', r_latent_specific, ...
    'dissipation_specific_Jkg', dissipation_specific, ...
    'strain_work_MPa', strain_work_MPa, ...
    'heat_source_Wkg', (latent_specific + dissipation_specific) / dt, ...
    'dissipation_nonnegative', dissipation_specific >= -1.0e-12);
end

function S = numerical_compliance(p, state, tau, T, dt)
base = sma.tensor_to_voigt(tau);
S = zeros(6);
for k = 1:6
    step = max(1.0e-4, 1.0e-6 * max(abs(base(k)), 1.0));
    plus = base;
    minus = base;
    plus(k) = plus(k) + step;
    minus(k) = minus(k) - step;
    [sp, ~] = sma.update_3d_stress( ...
        p, state, sma.voigt_to_tensor(plus), T, dt, false);
    [sm, ~] = sma.update_3d_stress( ...
        p, state, sma.voigt_to_tensor(minus), T, dt, false);
    S(:, k) = (sma.tensor_to_voigt(sp.h) - ...
        sma.tensor_to_voigt(sm.h)) / (2.0 * step);
end
end
