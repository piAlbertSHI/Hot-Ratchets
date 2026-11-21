function [trial, aux] = update_1d_stress(p, state, sigma, T, dt, compute_tangent)
%SMA.UPDATE_1D_STRESS Implicit stress-controlled constitutive update.
%
% STATE is never modified. The returned TRIAL state is committed only by
% the caller after its outer equilibrium iteration converges.

if nargin < 6
    compute_tangent = true;
end
if nargin < 5 || isempty(dt)
    dt = p.default_dt;
end
dt = max(dt, eps);
assert(isfinite(sigma) && isfinite(T), 'Stress and temperature must be finite.');
state = canonicalize_derived_state(p, state);
state=sma.ensure_zhang_direct_state(p,state);
state=sma.ensure_r_phase_state(p,state);

[phi_f, phi_r] = sma.trial_surfaces_1d(p,state,sigma,T);
branch = choose_branch(p, state, sigma, T, phi_f, phi_r);

if branch == 0
    trial = elastic_state(p, state, sigma, T);
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
    [x0, ~] = initial_guess(p, state, sigma, T, branch);
    if p.enable_composition_weighted_mt
        x0(1)=active_surface_fraction_seed( ...
            p,state,sigma,T,branch,reverse,x0(1));
    end
    boundary_probe = x0;
    boundary_probe(1) = double(branch > 0);
    boundary_R = sma.residual_1d( ...
        p,state,sigma,T,branch,reverse,boundary_probe);
    boundary = boundary_R(1) > p.numerics.surface_tol;
    if boundary
        [x, converged, iterations, residual_norm] = solve_boundary( ...
            p, state, sigma, T, branch, reverse, x0);
    else
        [x, converged, iterations, residual_norm] = solve_active( ...
            p, state, sigma, T, branch, reverse, x0);
        if ~converged
            [boundary_x, boundary_converged, boundary_iterations, ...
                boundary_residual] = ...
                solve_boundary(p, state, sigma, T, branch, reverse, x);
            boundary_R = sma.residual_1d( ...
                p,state,sigma,T,branch,reverse,boundary_x);
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
        error('sma:update_1d_stress:LocalNewton', ...
            'Local Newton failed: branch=%d, residual=%g.', branch, residual_norm);
    end
    [~, data] = sma.residual_1d(p, state, sigma, T, branch, reverse, x);
    trial = state;
    trial.T = T;
    trial.xi = min(max(x(1), 0.0), 1.0);
    trial.eps_tr = x(2);
    trial.eps_tp = x(3);
    trial.zeta_d = max(x(4), state.zeta_d);
    trial.beta = x(5);
    trial.sigma = sigma;
    trial.branch = branch;
    trial.direction = data.direction;
    trial.xi_d = data.xi_d;
    trial.Hcur = data.Hcur;
    trial.Hmax = data.Hmax;
    trial.xi_res = state.xi_res;
    if branch > 0
        trial.xi_r = trial.xi;
        trial.eps_tr_r = trial.eps_tr;
    else
        trial.xi_res = sma.residual_martensite_limit(p,state);
        if trial.xi <= trial.xi_res+p.numerics.xi_tol
            trial.xi = trial.xi_res;
        end
    end
    trial.Y_population_activity=sma.update_population_activity( ...
        p,state.Y_population_activity,state.xi,trial.xi);
    history_final=sma.candidate_zhang_history_1d( ...
        p,state,trial.xi,trial.eps_tp,sigma);
    trial=sma.commit_zhang_direct_history(trial,history_final);
    trial.eps_arp = sma.update_arp_1d(p,state,sigma,trial.xi);
    trial.zeta_p = sma.update_training_activity(p,state,trial,branch);
    E = sma.effective_modulus(p,trial.xi);
    alpha = sma.thermal_expansion(p,trial.xi);
    [eps_R,xi_R_trial] = sma.r_phase_strain_1d( ...
        p,sigma,T,trial.zeta_p,trial.xi,state.xi_R,state.xi);
    if p.enable_explicit_r_phase_state, trial.xi_R=xi_R_trial; end
    trial.eps = sigma/E+alpha*(T-p.T0)+eps_R+ ...
        trial.eps_tr+trial.eps_tp+trial.eps_arp;
    thermodynamic_target_xi = trial.xi;
    trial_surface_overdrive_MPa = choose(branch > 0, phi_f, phi_r);
    [trial, kinetic_rate, kinetic_fraction, activation_factor] = ...
        apply_rate_regularization( ...
        p,state,trial,branch,sigma,T,dt,trial_surface_overdrive_MPa);
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
    dsigma = max(1.0e-4, 1.0e-6 * max(abs(sigma), 1.0));
    [plus, ~] = sma.update_1d_stress(p, state, sigma + dsigma, T, dt, false);
    [minus, ~] = sma.update_1d_stress(p, state, sigma - dsigma, T, dt, false);
    aux.deps_dsigma = (plus.eps - minus.eps) / (2.0 * dsigma);
else
    aux.deps_dsigma = NaN;
end
end

function [trial, rate_1_per_s, kinetic_fraction, activation_factor] = ...
        apply_rate_regularization( ...
        p,state,trial,branch,sigma,T,dt,trial_surface_overdrive_MPa)
rate_1_per_s = inf;
kinetic_fraction = 1.0;
activation_factor = 1.0;
if ~p.enable_rate_regularization
    return;
end
if branch>0
    base_rate = p.phase_loading_rate_1_per_s;
else
    base_rate = p.phase_unloading_rate_1_per_s;
end
assert(isfinite(base_rate) && base_rate>0.0, ...
    'sma:update_1d_stress:InvalidRateRegularization', ...
    'The selected finite directional rate must be positive and finite.');
activation_factor = sma.phase_activation_factor( ...
    p, trial_surface_overdrive_MPa);
assert(isfield(p, 'phase_rate_multiplier_mode') && strcmp( ...
    char(string(p.phase_rate_multiplier_mode)), ...
    'none_qian_eq3_direction_coefficient_exact'), ...
    'sma:update_1d_stress:UnsupportedPhaseRateMultiplier', ...
    'The unique Qian Eq. (3) model requires an exact directional rate.');
% Qian Eq. (3) uses a directional material coefficient. The Gibbs return
% map already selects the branch and target, so the surface-distance
% diagnostic must not multiply that coefficient a second time.
rate_1_per_s = base_rate;
kinetic_fraction = -expm1(-rate_1_per_s*dt);
kinetic_fraction = min(max(kinetic_fraction,0.0),1.0);
if kinetic_fraction >= 1.0-10.0*eps
    return;
end

target = trial;
trial.xi = state.xi + kinetic_fraction*(target.xi-state.xi);
trial.eps_tr = state.eps_tr + ...
    kinetic_fraction*(target.eps_tr-state.eps_tr);
trial.eps_tp = state.eps_tp + ...
    kinetic_fraction*(target.eps_tp-state.eps_tp);
trial.zeta_d = max(state.zeta_d,state.zeta_d + ...
    kinetic_fraction*(target.zeta_d-state.zeta_d));
trial.beta = state.beta + kinetic_fraction*(target.beta-state.beta);
trial.xi_res = target.xi_res;
if branch>0
    trial.xi_r = trial.xi;
    trial.eps_tr_r = trial.eps_tr;
else
    trial.xi_r = state.xi_r;
    trial.eps_tr_r = state.eps_tr_r;
end
trial.Y_population_activity=sma.update_population_activity( ...
    p,state.Y_population_activity,state.xi,trial.xi);
history_final=sma.candidate_zhang_history_1d( ...
    p,state,trial.xi,trial.eps_tp,sigma);
trial=sma.commit_zhang_direct_history(trial,history_final);
trial.eps_arp = sma.update_arp_1d(p,state,sigma,trial.xi);
trial.zeta_p = sma.update_training_activity(p,state,trial,branch);
[trial.Hcur,trial.Hmax] = sma.hcur( ...
    p,abs(sigma+trial.beta),trial.zeta_d);
trial.xi_d = min(max(trial.Hcur/max(trial.Hmax,eps),0.0),1.0)*trial.xi;
E = sma.effective_modulus(p,trial.xi);
alpha = sma.thermal_expansion(p,trial.xi);
[eps_R,xi_R_trial] = sma.r_phase_strain_1d( ...
    p,sigma,T,trial.zeta_p,trial.xi,state.xi_R,state.xi);
if p.enable_explicit_r_phase_state, trial.xi_R=xi_R_trial; end
trial.eps = sigma/E+alpha*(T-p.T0)+eps_R+ ...
    trial.eps_tr+trial.eps_tp+trial.eps_arp;
end

function branch = choose_branch(p, state, sigma, T, phi_f, phi_r)
tol = p.numerics.surface_tol;
loading_measure = abs(sigma + state.beta) - abs(state.sigma + state.beta);
temperature_change = T - state.T;
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
    elseif can_reverse
        branch = -1;
    else
        branch = 0;
    end
elseif loading_measure < -mechanical_tol
    if can_reverse
        branch = -1;
    elseif can_forward
        branch = 1;
    else
        branch = 0;
    end
else
    if can_forward && can_reverse
        branch = choose(temperature_change < 0,1,-1);
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
reverse = struct('xi_r', state.xi_r, 'eps_tr_r', state.eps_tr_r);
if branch < 0 && reverse.xi_r <= 0.0
    reverse.xi_r = state.xi;
    reverse.eps_tr_r = state.eps_tr;
end
end

function [x0, boundary] = initial_guess(p, state, sigma, T, branch)
sigma_eff = abs(sigma + state.beta);
if branch > 0
    xi_target = (sigma_eff / p.CM - (T - p.Ms)) / (p.Ms - p.Mf);
    xi_target = max(xi_target, state.xi);
    boundary = xi_target >= 1.0;
else
    xi_target = (sigma_eff / p.CA - (T - p.Af)) / (p.Af - p.As);
    xi_target = min(xi_target, state.xi);
    reverse_limit = sma.residual_martensite_limit(p,state);
    boundary = xi_target <= reverse_limit;
    xi_target = max(xi_target,reverse_limit);
end
xi_target = min(max(xi_target, 0.0), 1.0);
x0 = [xi_target; state.eps_tr; state.eps_tp; state.zeta_d; state.beta];
end

function xi_seed=active_surface_fraction_seed( ...
    p,state,sigma,T,branch,reverse,xi_fallback)
% Bracket only the constitutive surface to avoid a singular Newton start at
% the ninth-order hardening origin. Memory variables remain committed here
% and are solved simultaneously by solve_active afterwards.
if branch>0
    lo=state.xi; hi=1.0;
else
    lo=sma.residual_martensite_limit(p,state); hi=state.xi;
end
surface=@(xi)surface_at_fraction( ...
    p,state,sigma,T,branch,reverse,xi);
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

function value=surface_at_fraction( ...
    p,state,sigma,T,branch,reverse,xi)
x=[xi;state.eps_tr;state.eps_tp;state.zeta_d;state.beta];
R=sma.residual_1d(p,state,sigma,T,branch,reverse,x);
value=R(1);
end

function [x, converged, iterations, residual_norm] = solve_active( ...
    p, state, sigma, T, branch, reverse, x)
y = [x(1);state.zeta_d;state.beta];
fun = @(memory) active_memory_residual( ...
    p,state,sigma,T,branch,reverse,memory);
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
    lower = [xi_lower;state.zeta_d;-inf];
    upper = [xi_upper;inf;inf];
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
x_probe = [y(1);state.eps_tr;state.eps_tp;y(2);y(3)];
[~,data] = sma.residual_1d(p,state,sigma,T,branch,reverse,x_probe);
delta_xi = y(1)-state.xi;
x = [y(1);state.eps_tr+data.lambda_tr*delta_xi; ...
    state.eps_tp+data.lambda_tp*delta_xi;y(2);y(3)];
final_R = sma.residual_1d(p,state,sigma,T,branch,reverse,x);
residual_norm = norm(final_R,inf);
converged = converged || residual_norm<=p.numerics.residual_tol;
end

function [x, converged, iterations, residual_norm] = solve_boundary( ...
    p, state, sigma, T, branch, reverse, ~)
if branch>0
    fixed_xi = 1.0;
else
    fixed_xi = sma.residual_martensite_limit(p,state);
end
y = [state.zeta_d;state.beta];
converged = false;
for iterations = 1:100
    x_probe = [fixed_xi;state.eps_tr;state.eps_tp;y(1);y(2)];
    [~,data] = sma.residual_1d(p,state,sigma,T,branch,reverse,x_probe);
    target = [state.zeta_d+abs(data.xi_d-state.xi_d);data.beta_target];
    residual_norm = norm(y-target,inf);
    if residual_norm<=p.numerics.residual_tol
        converged = true;
        break;
    end
    y = 0.5*y+0.5*target;
    y(1) = max(y(1),state.zeta_d);
end
x_probe = [fixed_xi;state.eps_tr;state.eps_tp;y(1);y(2)];
[~,data] = sma.residual_1d(p,state,sigma,T,branch,reverse,x_probe);
delta_xi = fixed_xi-state.xi;
x = [fixed_xi; ...
    state.eps_tr+data.lambda_tr*delta_xi; ...
    state.eps_tp+data.lambda_tp*delta_xi; ...
    y(1);y(2)];
final_R = sma.residual_1d(p,state,sigma,T,branch,reverse,x);
if branch<0
    residual_norm = norm(final_R(3:5),inf);
else
    residual_norm = norm(final_R(2:5),inf);
end
converged = converged || residual_norm<=p.numerics.residual_tol;
end

function R = active_memory_residual( ...
    p,state,sigma,T,branch,reverse,memory)
x = [memory(1);state.eps_tr;state.eps_tp;memory(2);memory(3)];
[all_R,data] = sma.residual_1d( ...
    p,state,sigma,T,branch,reverse,x);
R = all_R([1,4,5]);
delta_xi_d = data.xi_d-state.xi_d;
if branch > 0
    R(2) = memory(2)-state.zeta_d-delta_xi_d;
else
    R(2) = memory(2)-state.zeta_d+delta_xi_d;
end
end

function [candidate, accepted] = line_search_active(fun, x, step, state, branch, ...
    residual_norm, p)
eta = 1.0;
accepted = false;
candidate = x;
for ls = 1:p.numerics.max_line_search
    value = x + eta * step;
    value(1) = min(max(value(1),0.0),1.0);
    value(2) = max(value(2),state.zeta_d);
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

function trial = elastic_state(p, state, sigma, T)
trial = state;
trial.T = T;
trial.sigma = sigma;
trial.branch = 0.0;
[trial.Hcur,trial.Hmax] = sma.hcur( ...
    p,abs(sigma+state.beta),state.zeta_d);
trial.xi_d = min(max(trial.Hcur/max(trial.Hmax,eps),0.0),1.0) * ...
    state.xi;
trial.eps_arp = sma.update_arp_1d(p,state,sigma,state.xi);
E = sma.effective_modulus(p, state.xi);
alpha = sma.thermal_expansion(p, state.xi);
[eps_R,xi_R_trial] = sma.r_phase_strain_1d( ...
    p,sigma,T,state.zeta_p,state.xi,state.xi_R,state.xi);
if p.enable_explicit_r_phase_state, trial.xi_R=xi_R_trial; end
trial.eps = sigma / E + alpha * (T - p.T0) + ...
    eps_R + state.eps_tr + state.eps_tp + trial.eps_arp;
end

function state = canonicalize_derived_state(p, state)
state=sma.ensure_population_state(p,state);
state=sma.ensure_r_phase_state(p,state);
[state.Hcur,state.Hmax] = sma.hcur( ...
    p,abs(state.sigma+state.beta),state.zeta_d);
state.xi_d = min(max(state.Hcur/max(state.Hmax,eps),0.0),1.0) * ...
    state.xi;
end

function aux = make_aux(p, state, trial, dt, converged, iterations, ...
    residual_norm, phi)
delta_xi = trial.xi - state.xi;
delta_eps_tr = trial.eps_tr - state.eps_tr;
delta_eps_tp = trial.eps_tp - state.eps_tp;
delta_eps_arp = trial.eps_arp - state.eps_arp;
latent_value = sma.effective_latent_heat(p,state,trial);
latent_specific = latent_value * delta_xi;
r_latent_specific = sma.r_phase_latent_specific( ...
    p,state.xi,abs(state.sigma),state.T,state.zeta_p, ...
    trial.xi,abs(trial.sigma),trial.T,trial.zeta_p, ...
    state.xi_R,trial.xi_R);
latent_specific = latent_specific+r_latent_specific;
mechanical_specific = 1.0e6 / p.rho * abs( ...
    (trial.sigma + trial.beta) * delta_eps_tr + ...
    trial.sigma * (delta_eps_tp+delta_eps_arp));
dissipation_specific = p.taylor_quinney * mechanical_specific;
strain_work_MPa = 0.5*(state.sigma+trial.sigma)* ...
    (trial.eps-state.eps);
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
