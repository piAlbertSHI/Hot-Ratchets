function [trial,aux] = monolithic_coupled_step_1d( ...
    p,state,mode,target,Tamb,dt,heat_loss_rate)
%SMA.MONOLITHIC_COUPLED_STEP_1D Joint phase/temperature return mapping.
% The routine is used only after the nested thermo-mechanical fixed point
% encounters an active-set cycle.  Temperature, phase fraction and local
% memory are solved in one residual system; no two branch states are averaged.

assert(mode==0 || mode==1,'Monolithic coupling supports stress/strain control.');
state=sma.ensure_population_state(p,state);
state=sma.ensure_zhang_direct_state(p,state);
state=sma.ensure_r_phase_state(p,state);
if mode==1
    measure=target-state.eps;
else
    measure=target-state.sigma;
end
if measure>=-p.numerics.direction_tol
    branches=[1,-1];
else
    branches=[-1,1];
end
best_norm=inf;
best=[];
for branch=branches
    lower_xi=sma.residual_martensite_limit(p,state);
    if branch>0
        if state.xi>=1.0-p.numerics.xi_tol, continue; end
        xi0=min(state.xi+max(1.0e-5,0.01*(1.0-state.xi)),1.0);
    else
        if state.xi<=lower_xi+p.numerics.xi_tol, continue; end
        xi0=max(state.xi-max(1.0e-5,0.01*(state.xi-lower_xi)),lower_xi);
    end
    temperature0=state.T-heat_loss_rate*(state.T-Tamb)*dt/p.cp;
    if mode==1
        sigma0=state.sigma+sma.effective_modulus(p,state.xi)*measure;
        y=[sigma0;xi0;state.zeta_d;state.beta;temperature0];
    else
        y=[xi0;state.zeta_d;state.beta;temperature0];
    end
    fun=@(value)evaluate_residual( ...
        p,state,mode,target,Tamb,dt,heat_loss_rate,branch,value);
    [y,converged,iterations,residual_norm]=solve_newton( ...
        fun,y,p,state,branch,mode);
    if residual_norm<best_norm
        best_norm=residual_norm;
        best=struct('y',y,'branch',branch,'iterations',iterations, ...
            'converged',converged);
    end
    if converged, break; end
end
if isempty(best) || ~best.converged
    error('sma:monolithic_coupled_step_1d:Newton', ...
        'Joint thermo-mechanical return mapping failed: residual %.6g.',best_norm);
end
[~,trial,aux]=evaluate_residual( ...
    p,state,mode,target,Tamb,dt,heat_loss_rate,best.branch,best.y);
aux.thermal_iterations=best.iterations;
aux.thermal_converged=true;
aux.thermal_branch_switches=1;
aux.thermal_cutbacks=0;
aux.monolithic_thermal_solves=1;
end

function [R,trial,aux]=evaluate_residual( ...
    p,state,mode,target,Tamb,dt,heat_loss_rate,branch,y)
if mode==1
    sigma=y(1); memory=y(2:4); T=y(5);
else
    sigma=target; memory=y(1:3); T=y(4);
end
xi=memory(1); zeta_d=memory(2); beta=memory(3);
reverse=reverse_reference(state,branch);
x_probe=[xi;state.eps_tr;state.eps_tp;zeta_d;beta];
[~,data]=sma.residual_1d(p,state,sigma,T,branch,reverse,x_probe);
delta_xi=xi-state.xi;
x=[xi;state.eps_tr+data.lambda_tr*delta_xi; ...
    state.eps_tp+data.lambda_tp*delta_xi;zeta_d;beta];
[all_R,data]=sma.residual_1d(p,state,sigma,T,branch,reverse,x);

trial=state;
trial.T=T;
trial.xi=min(max(x(1),0.0),1.0);
trial.eps_tr=x(2);
trial.eps_tp=x(3);
trial.zeta_d=max(x(4),state.zeta_d);
trial.beta=x(5);
trial.sigma=sigma;
trial.branch=branch;
trial.direction=data.direction;
trial.xi_d=data.xi_d;
trial.Hcur=data.Hcur;
trial.Hmax=data.Hmax;
trial.xi_res=state.xi_res;
if branch>0
    trial.xi_r=trial.xi;
    trial.eps_tr_r=trial.eps_tr;
else
    trial.xi_res=sma.residual_martensite_limit(p,state);
end
trial.Y_population_activity=sma.update_population_activity( ...
    p,state.Y_population_activity,state.xi,trial.xi);
history=sma.candidate_zhang_history_1d( ...
    p,state,trial.xi,trial.eps_tp,trial.sigma);
trial=sma.commit_zhang_direct_history(trial,history);
trial.zeta_p=sma.update_training_activity(p,state,trial,branch);
trial.eps_arp=sma.update_arp_1d(p,state,trial.sigma,trial.xi);
E=sma.effective_modulus(p,trial.xi);
alpha=sma.thermal_expansion(p,trial.xi);
[eps_R,xi_R]=sma.r_phase_strain_1d( ...
    p,trial.sigma,T,trial.zeta_p,trial.xi,state.xi_R,state.xi);
if p.enable_explicit_r_phase_state, trial.xi_R=xi_R; end
trial.eps=trial.sigma/E+alpha*(T-p.T0)+eps_R+ ...
    trial.eps_tr+trial.eps_tp+trial.eps_arp;

aux=make_aux(p,state,trial,dt,data.phi);
heat_loss=heat_loss_rate*(T-Tamb)*dt;
thermal_residual_Jkg=p.cp*(T-state.T)- ...
    (aux.latent_specific_Jkg+aux.dissipation_specific_Jkg-heat_loss);
thermal_residual_K=thermal_residual_Jkg/p.cp;
aux.heat_loss_specific_Jkg=heat_loss;
aux.heat_loss_rate_WkgK=heat_loss_rate;
aux.temperature_change=T-state.T;
aux.thermomechanical_temperature_residual_K=thermal_residual_K;
if mode==1
    control_residual=(trial.eps-target)/max(p.Hmax_i,1.0e-4);
    R=[control_residual;all_R([1,4,5]);thermal_residual_K];
else
    R=[all_R([1,4,5]);thermal_residual_K];
end
aux.residual_norm=norm(all_R,inf);
end

function [y,converged,iterations,residual_norm]= ...
    solve_newton(fun,y,p,state,branch,mode)
converged=false;
for iterations=1:2*p.numerics.max_iterations
    R=fun(y);
    residual_norm=norm(R,inf);
    if residuals_converged(R,p)
        converged=true;
        break;
    end
    J=sma.numerical_jacobian(fun,y,p.numerics.jacobian_step);
    step=-J\R;
    [candidate,accepted]=line_search( ...
        fun,y,step,residual_norm,p,state,branch,mode);
    if ~accepted
        step=-(J.'*J+1.0e-8*eye(size(J,2)))\(J.'*R);
        [candidate,accepted]=line_search( ...
            fun,y,step,residual_norm,p,state,branch,mode);
    end
    if ~accepted || norm(candidate-y,inf)<=p.numerics.step_tol
        break;
    end
    y=candidate;
end
residual_norm=norm(fun(y),inf);
converged=converged || residuals_converged(fun(y),p);
end

function pass=residuals_converged(R,p)
pass=norm(R(1:end-1),inf)<=p.numerics.residual_tol && ...
    abs(R(end))<=p.thermomechanical_coupling_temperature_tol_K;
end

function [candidate,accepted]=line_search( ...
    fun,y,step,residual_norm,p,state,branch,mode)
candidate=y;
accepted=false;
for k=1:2*p.numerics.max_line_search
    value=y+(0.5^(k-1))*step;
    if mode==1
        value(1)=min(max(value(1),-p.numerics.stress_solve_bound), ...
            p.numerics.stress_solve_bound);
        xi_index=2; zeta_index=3; temperature_index=5;
    else
        xi_index=1; zeta_index=2; temperature_index=4;
    end
    lower_xi=sma.residual_martensite_limit(p,state);
    value(xi_index)=min(max(value(xi_index),lower_xi),1.0);
    value(zeta_index)=max(value(zeta_index),state.zeta_d);
    value(temperature_index)=min(max(value(temperature_index), ...
        p.min_temperature),p.max_temperature);
    if branch>0
        value(xi_index)=max(value(xi_index),state.xi);
    else
        value(xi_index)=min(value(xi_index),state.xi);
    end
    if norm(fun(value),inf)<residual_norm
        candidate=value;
        accepted=true;
        return;
    end
end
end

function reverse=reverse_reference(state,branch)
reverse=struct('xi_r',state.xi_r,'eps_tr_r',state.eps_tr_r);
if branch<0 && reverse.xi_r<=0.0
    reverse.xi_r=state.xi;
    reverse.eps_tr_r=state.eps_tr;
end
end

function aux=make_aux(p,state,trial,dt,phi)
delta_xi=trial.xi-state.xi;
delta_eps_tr=trial.eps_tr-state.eps_tr;
delta_eps_tp=trial.eps_tp-state.eps_tp;
delta_eps_arp=trial.eps_arp-state.eps_arp;
latent_value=sma.effective_latent_heat(p,state,trial);
latent_specific=latent_value*delta_xi;
r_latent_specific=sma.r_phase_latent_specific( ...
    p,state.xi,abs(state.sigma),state.T,state.zeta_p, ...
    trial.xi,abs(trial.sigma),trial.T,trial.zeta_p, ...
    state.xi_R,trial.xi_R);
latent_specific=latent_specific+r_latent_specific;
mechanical_specific=1.0e6/p.rho*abs( ...
    (trial.sigma+trial.beta)*delta_eps_tr+ ...
    trial.sigma*(delta_eps_tp+delta_eps_arp));
dissipation_specific=p.taylor_quinney*mechanical_specific;
aux=struct('converged',true,'iterations',0,'residual_norm',inf, ...
    'phi',phi,'branch',trial.branch,'delta_xi',delta_xi, ...
    'effective_latent_Jkg',latent_value, ...
    'latent_specific_Jkg',latent_specific, ...
    'r_phase_latent_specific_Jkg',r_latent_specific, ...
    'dissipation_specific_Jkg',dissipation_specific, ...
    'strain_work_MPa',0.5*(state.sigma+trial.sigma)*(trial.eps-state.eps), ...
    'heat_source_Wkg',(latent_specific+dissipation_specific)/max(dt,eps), ...
    'dissipation_nonnegative',dissipation_specific>=-1.0e-12);
end
