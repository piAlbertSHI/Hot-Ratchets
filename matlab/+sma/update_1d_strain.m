function [trial, aux] = update_1d_strain(p, state, eps_target, T, dt, compute_tangent)
%SMA.UPDATE_1D_STRAIN Outer stress solve with rollback-safe trial updates.

if nargin < 6
    compute_tangent = true;
end
if nargin < 5 || isempty(dt)
    dt = p.default_dt;
end
state=sma.ensure_population_state(p,state);
state=sma.ensure_zhang_direct_state(p,state);
state=sma.ensure_r_phase_state(p,state);
bound = p.numerics.stress_solve_bound;
residual = @(sigma) strain_residual(p, state, sigma, eps_target, T, dt);
lo = -bound;
hi = bound;
rlo = residual(lo);
rhi = residual(hi);
if rlo * rhi <= 0.0
    sigma = fzero(residual, [lo, hi]);
else
    sigma = fzero(residual, state.sigma);
end
try
    [trial, aux] = sma.update_1d_stress( ...
        p, state, sigma, T, dt, compute_tangent);
catch ME
    if p.numerics.cutback_on_failure && strcmp(ME.identifier, ...
            'sma:update_1d_stress:LocalNewton')
        [trial,aux] = solve_direct_strain( ...
            p,state,eps_target,T,dt,sigma);
        return;
    end
    rethrow(ME);
end
aux.strain_residual = trial.eps - eps_target;
strain_tolerance = max(1.0e-9, ...
    p.numerics.residual_tol * max(p.Hmax_i, 1.0e-4));
if abs(aux.strain_residual) > strain_tolerance
    if p.numerics.cutback_on_failure
        [trial, aux] = solve_direct_strain( ...
            p,state,eps_target,T,dt,sigma);
    else
        try
            [trial, aux] = solve_direct_strain( ...
                p,state,eps_target,T,dt,sigma);
        catch ME
            if ~strcmp(ME.identifier, ...
                    'sma:update_1d_strain:CoupledNewton')
                rethrow(ME);
            end
        end
    end
end
end

function value = strain_residual(p, state, sigma, eps_target, T, dt)
trial_p = p;
% Root bracketing intentionally probes remote stresses. Those trial points
% are not committed and may lie outside an admissible active surface.
trial_p.numerics.cutback_on_failure = false;
[trial, ~] = sma.update_1d_stress(trial_p, state, sigma, T, dt, false);
value = trial.eps - eps_target;
end

function [trial, aux] = solve_direct_strain( ...
    p,state,eps_target,T,dt,sigma_guess)
if eps_target >= state.eps
    branches = [1,-1];
else
    branches = [-1,1];
end
best = struct('converged',false,'norm',inf,'y',[],'branch',0, ...
    'iterations',0,'reverse',struct());
for branch = branches
    if (branch>0 && state.xi>=1.0-p.numerics.xi_tol) || ...
            (branch<0 && state.xi<= ...
            sma.residual_martensite_limit(p,state)+p.numerics.xi_tol)
        continue;
    end
    reverse = reverse_reference(state,branch);
    y = initial_direct_guess( ...
        p,state,eps_target,T,sigma_guess,branch);
    fun = @(value) direct_residual( ...
        p,state,eps_target,T,branch,reverse,value);
    [y,converged,iterations,residual_norm] = ...
        direct_newton(fun,y,state,branch,p);
    if residual_norm<best.norm
        best = struct('converged',converged,'norm',residual_norm, ...
            'y',y,'branch',branch,'iterations',iterations, ...
            'reverse',reverse);
    end
    if converged
        break;
    end
end
if ~best.converged
    for branch = branches
        reverse = reverse_reference(state,branch);
        [boundary_y,converged,iterations,residual_norm] = ...
            solve_direct_boundary( ...
            p,state,eps_target,T,branch,reverse,sigma_guess);
        if branch>0
            fixed_xi = 1.0;
        else
            fixed_xi = sma.residual_martensite_limit(p,state);
        end
        full_y = [boundary_y(1);fixed_xi;boundary_y(2:3)];
        [boundary_x,~] = reconstruct_state( ...
            p,state,full_y(1),T,branch,reverse,full_y(2:4));
        boundary_R = sma.residual_1d( ...
            p,state,full_y(1),T,branch,reverse,boundary_x);
        admissible = boundary_R(1)>=-p.numerics.surface_tol;
        if converged && admissible
            best = struct('converged',true,'norm',residual_norm, ...
                'y',full_y,'branch',branch,'iterations',iterations, ...
                'reverse',reverse);
            break;
        end
        if residual_norm<best.norm
            best.norm = residual_norm;
        end
    end
end
if ~best.converged
    error('sma:update_1d_strain:CoupledNewton', ...
        ['Coupled strain return mapping failed: residual=%g, ', ...
        'branch=%d, sigma=%g, xi=%g, zeta=%g, beta=%g.'], ...
        best.norm,best.branch,best.y(1),best.y(2),best.y(3),best.y(4));
end

function [y,converged,iterations,residual_norm] = ...
    solve_direct_boundary( ...
    p,state,eps_target,T,branch,reverse,sigma_guess)
y = [sigma_guess;state.zeta_d;state.beta];
fun = @(value) boundary_strain_residual( ...
    p,state,eps_target,T,branch,reverse,value);
converged = false;
for iterations = 1:2*p.numerics.max_iterations
    R = fun(y);
    residual_norm = norm(R,inf);
    if residual_norm<=p.numerics.residual_tol
        converged = true;
        break;
    end
    J = sma.numerical_jacobian(fun,y,p.numerics.jacobian_step);
    step = -J\R;
    eta = 1.0;
    accepted = false;
    for idx = 1:2*p.numerics.max_line_search
        candidate = y+eta*step;
        candidate(1) = min(max(candidate(1), ...
            -p.numerics.stress_solve_bound),p.numerics.stress_solve_bound);
        candidate(2) = max(candidate(2),state.zeta_d);
        if norm(fun(candidate),inf)<residual_norm
            y = candidate;
            accepted = true;
            break;
        end
        eta = 0.5*eta;
    end
    if ~accepted
        break;
    end
end
residual_norm = norm(fun(y),inf);
converged = converged || residual_norm<=p.numerics.residual_tol;
end

function R = boundary_strain_residual( ...
    p,state,eps_target,T,branch,reverse,y)
if branch>0
    fixed_xi = 1.0;
else
    fixed_xi = sma.residual_martensite_limit(p,state);
end
memory = [fixed_xi;y(2);y(3)];
[x,~] = reconstruct_state(p,state,y(1),T,branch,reverse,memory);
all_R = sma.residual_1d(p,state,y(1),T,branch,reverse,x);
E = sma.effective_modulus(p,x(1));
alpha = sma.thermal_expansion(p,x(1));
eps_R = sma.r_phase_strain_1d( ...
    p,y(1),T,state.zeta_p,x(1),state.xi_R,state.xi);
eps_arp = sma.update_arp_1d(p,state,y(1),x(1));
eps_value = y(1)/E+alpha*(T-p.T0)+eps_R+x(2)+x(3)+eps_arp;
if branch<0
    memory_residual = all_R([3,4,5]);
    memory_residual = memory_residual(2:3);
else
    memory_residual = all_R([4,5]);
end
R = [(eps_value-eps_target)/max(p.Hmax_i,1.0e-4); ...
    memory_residual];
end
[x,data] = reconstruct_state( ...
    p,state,best.y(1),T,best.branch,best.reverse,best.y(2:4));
trial = state;
trial.T = T;
trial.xi = min(max(x(1),0.0),1.0);
trial.eps_tr = x(2);
trial.eps_tp = x(3);
trial.zeta_d = max(x(4),state.zeta_d);
trial.beta = x(5);
trial.sigma = best.y(1);
trial.branch = best.branch;
trial.direction = data.direction;
trial.xi_d = data.xi_d;
trial.Hcur = data.Hcur;
trial.Hmax = data.Hmax;
trial.xi_res = state.xi_res;
if best.branch>0
    trial.xi_r = trial.xi;
    trial.eps_tr_r = trial.eps_tr;
else
    trial.xi_res = sma.residual_martensite_limit(p,state);
    if trial.xi<=trial.xi_res+p.numerics.xi_tol
        trial.xi = trial.xi_res;
    end
end
trial.Y_population_activity=sma.update_population_activity( ...
    p,state.Y_population_activity,state.xi,trial.xi);
history_final=sma.candidate_zhang_history_1d( ...
    p,state,trial.xi,trial.eps_tp,trial.sigma);
trial=sma.commit_zhang_direct_history(trial,history_final);
trial.zeta_p = sma.update_training_activity( ...
    p,state,trial,best.branch);
trial.eps_arp = sma.update_arp_1d(p,state,trial.sigma,trial.xi);
E = sma.effective_modulus(p,trial.xi);
alpha = sma.thermal_expansion(p,trial.xi);
[eps_R,xi_R_trial] = sma.r_phase_strain_1d( ...
    p,trial.sigma,T,trial.zeta_p,trial.xi,state.xi_R,state.xi);
if p.enable_explicit_r_phase_state, trial.xi_R=xi_R_trial; end
trial.eps = trial.sigma/E+alpha*(T-p.T0)+ ...
    eps_R+trial.eps_tr+trial.eps_tp+trial.eps_arp;
aux = direct_aux(p,state,trial,dt,best.iterations,best.norm,data.phi);
aux.strain_residual = trial.eps-eps_target;
aux.deps_dsigma = NaN;
end

function reverse = reverse_reference(state,branch)
reverse = struct('xi_r',state.xi_r,'eps_tr_r',state.eps_tr_r);
if branch<0 && reverse.xi_r<=0.0
    reverse.xi_r = state.xi;
    reverse.eps_tr_r = state.eps_tr;
end
end

function y = initial_direct_guess( ...
    p,state,eps_target,T,sigma_guess,branch)
E = sma.effective_modulus(p,state.xi);
alpha = sma.thermal_expansion(p,state.xi);
elastic_memory_eps = sigma_guess/E+alpha*(T-p.T0)+ ...
    state.eps_tr+state.eps_tp+state.eps_arp;
[Hcur,~] = sma.hcur(p,abs(sigma_guess+state.beta),state.zeta_d);
delta_xi = (eps_target-elastic_memory_eps)/max(Hcur,1.0e-4);
if branch>0
    xi = max(state.xi,state.xi+delta_xi);
else
    xi = min(state.xi,state.xi+delta_xi);
    xi = max(xi,sma.residual_martensite_limit(p,state));
end
xi = min(max(xi,0.0),1.0);
y = [sigma_guess;xi;state.zeta_d;state.beta];
end

function [y,converged,iterations,residual_norm] = ...
    direct_newton(fun,y,state,branch,p)
converged = false;
for iterations = 1:2*p.numerics.max_iterations
    R = fun(y);
    residual_norm = norm(R,inf);
    if residual_norm<=p.numerics.residual_tol
        converged = true;
        break;
    end
    J = sma.numerical_jacobian(fun,y,p.numerics.jacobian_step);
    step = -J\R;
    [candidate,accepted] = direct_line_search( ...
        fun,y,step,state,branch,residual_norm,p);
    if ~accepted
        step = -(J.'*J+1.0e-8*eye(4))\(J.'*R);
        [candidate,accepted] = direct_line_search( ...
            fun,y,step,state,branch,residual_norm,p);
    end
    if ~accepted || norm(candidate-y,inf)<=p.numerics.step_tol
        break;
    end
    y = candidate;
end
residual_norm = norm(fun(y),inf);
converged = converged || residual_norm<=p.numerics.residual_tol;
end

function [candidate,accepted] = direct_line_search( ...
    fun,y,step,state,branch,residual_norm,p)
accepted = false;
candidate = y;
eta = 1.0;
for idx = 1:2*p.numerics.max_line_search
    value = y+eta*step;
    value(1) = min(max(value(1), ...
        -p.numerics.stress_solve_bound),p.numerics.stress_solve_bound);
    value(2) = min(max(value(2),0.0),1.0);
    value(3) = max(value(3),state.zeta_d);
    if branch>0
        value(2) = max(value(2),state.xi);
    else
        reverse_limit = sma.residual_martensite_limit(p,state);
        value(2) = min(max(value(2),reverse_limit),state.xi);
    end
    if norm(fun(value),inf)<residual_norm
        candidate = value;
        accepted = true;
        return;
    end
    eta = 0.5*eta;
end
end

function R = direct_residual( ...
    p,state,eps_target,T,branch,reverse,y)
[x,~] = reconstruct_state(p,state,y(1),T,branch,reverse,y(2:4));
all_R = sma.residual_1d(p,state,y(1),T,branch,reverse,x);
E = sma.effective_modulus(p,x(1));
alpha = sma.thermal_expansion(p,x(1));
eps_R = sma.r_phase_strain_1d( ...
    p,y(1),T,state.zeta_p,x(1),state.xi_R,state.xi);
eps_arp = sma.update_arp_1d(p,state,y(1),x(1));
eps_value = y(1)/E+alpha*(T-p.T0)+eps_R+x(2)+x(3)+eps_arp;
R = [(eps_value-eps_target)/max(p.Hmax_i,1.0e-4); ...
    all_R([1,4,5])];
end

function [x,data] = reconstruct_state( ...
    p,state,sigma,T,branch,reverse,memory)
x_probe = [memory(1);state.eps_tr;state.eps_tp;memory(2);memory(3)];
[~,data] = sma.residual_1d( ...
    p,state,sigma,T,branch,reverse,x_probe);
delta_xi = memory(1)-state.xi;
x = [memory(1); ...
    state.eps_tr+data.lambda_tr*delta_xi; ...
    state.eps_tp+data.lambda_tp*delta_xi; ...
    memory(2);memory(3)];
end

function aux = direct_aux(p,state,trial,dt,iterations,residual_norm,phi)
delta_xi = trial.xi-state.xi;
delta_eps_tr = trial.eps_tr-state.eps_tr;
delta_eps_tp = trial.eps_tp-state.eps_tp;
delta_eps_arp = trial.eps_arp-state.eps_arp;
latent_value = sma.effective_latent_heat(p,state,trial);
latent_specific = latent_value*delta_xi;
r_latent_specific = sma.r_phase_latent_specific( ...
    p,state.xi,abs(state.sigma),state.T,state.zeta_p, ...
    trial.xi,abs(trial.sigma),trial.T,trial.zeta_p, ...
    state.xi_R,trial.xi_R);
latent_specific = latent_specific+r_latent_specific;
mechanical_specific = 1.0e6/p.rho*abs( ...
    (trial.sigma+trial.beta)*delta_eps_tr+ ...
    trial.sigma*(delta_eps_tp+delta_eps_arp));
dissipation_specific = p.taylor_quinney*mechanical_specific;
strain_work_MPa = 0.5*(state.sigma+trial.sigma)* ...
    (trial.eps-state.eps);
aux = struct('converged',true,'iterations',iterations, ...
    'residual_norm',residual_norm,'phi',phi,'branch',trial.branch, ...
    'delta_xi',delta_xi,'effective_latent_Jkg',latent_value, ...
    'latent_specific_Jkg',latent_specific, ...
    'r_phase_latent_specific_Jkg',r_latent_specific, ...
    'dissipation_specific_Jkg',dissipation_specific, ...
    'strain_work_MPa',strain_work_MPa, ...
    'heat_source_Wkg',(latent_specific+dissipation_specific)/max(dt,eps), ...
    'dissipation_nonnegative',dissipation_specific>=-1.0e-12);
end
