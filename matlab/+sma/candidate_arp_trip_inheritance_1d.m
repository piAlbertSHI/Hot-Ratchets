function [state,audit] = candidate_arp_trip_inheritance_1d( ...
    old,xi_M_new,sigma_MPa,dt,closure,guess)
%SMA.CANDIDATE_ARP_TRIP_INHERITANCE_1D Zhang ARP/TRIP shadow update.
% This source-faithful candidate implements Zhang (2022) Eqs. (8), (12)--
% (14), (23), (28), (40), (48), and (64)--(66).  It is deliberately
% isolated from the production return map and has no observation access.
%
% The phase change is applied first as an exactly conservative transfer of
% plastic strain, accumulated slip, and dislocation inventory.  Intrinsic
% ARP/TRIP is then integrated in the remaining A+R inventory.  The local
% slip laws are solved in the stable root form
%   pi_i - tau0*(Delta gamma_i/dt)^(1/mp) = 0,
% rather than evaluating the published 200th power directly.

required_state = {'xi_M','xi_M_acc','gamma_AR','gamma_M', ...
    'eps_p_AR','eps_p_M','rho_AR','rho_M','eps_ARP','eps_TRIP', ...
    'chi','beta','xi_Mr','Omega','Yp','DeltaEtaP_MPaK','Lambda_M'};
for k=1:numel(required_state)
    assert(isfield(old,required_state{k}), ...
        'sma:candidate_arp_trip_inheritance_1d:State', ...
        'Missing plastic state field %s.',required_state{k});
end
required_closure = {'chi_mu0','tau0_MPa','Hp_MPa','mp','chi_f0_MPa', ...
    'chi_p_sat_MPa','b','k1_per_m','k2','direction_MT', ...
    'allow_arp','allow_trip','max_local_slip_increment'};
for k=1:numel(required_closure)
    assert(isfield(closure,required_closure{k}), ...
        'sma:candidate_arp_trip_inheritance_1d:Closure', ...
        'Missing explicit closure field %s.',required_closure{k});
end
state_values = cellfun(@(name)old.(name),required_state);
closure_values = cellfun(@(name)double(closure.(name)),required_closure);
assert(all(isfinite([state_values,closure_values,xi_M_new,sigma_MPa,dt])) ...
    && dt>0, 'sma:candidate_arp_trip_inheritance_1d:Finite', ...
    'State, controls, time step, and closure values must be finite.');
tol=100*eps;
assert(old.xi_M>=-tol && old.xi_M<=1+tol && xi_M_new>=-tol && ...
    xi_M_new<=1+tol && old.xi_M_acc>=old.xi_M-tol && ...
    min([old.gamma_AR,old.gamma_M,old.rho_AR,old.rho_M])>=-tol, ...
    'sma:candidate_arp_trip_inheritance_1d:StateBounds', ...
    'Committed phase, slip, or dislocation state is inadmissible.');
assert(closure.chi_mu0>=0 && closure.tau0_MPa>0 && ...
    closure.Hp_MPa>0 && closure.mp>=1 && closure.chi_f0_MPa>=0 && ...
    closure.chi_p_sat_MPa>=0 && closure.b>=0 && ...
    closure.k1_per_m>=0 && closure.k2>=0 && ...
    any(closure.direction_MT==[-1,1]) && ...
    any(closure.allow_arp==[0,1]) && any(closure.allow_trip==[0,1]) && ...
    closure.max_local_slip_increment>0, ...
    'sma:candidate_arp_trip_inheritance_1d:ClosureRange', ...
    'ARP/TRIP closure values are inadmissible.');

xi_old=min(max(old.xi_M,0),1);
xi_new=min(max(xi_M_new,0),1);
ar_old=1-xi_old; ar_new=1-xi_new;
dxi=xi_new-xi_old;
dxi_plus=max(dxi,0); dxi_minus=max(-dxi,0);

% Conservative phase-interface transfer, before intrinsic plastic slip.
[P_M0,P_AR0]=transfer_inventory(xi_old,xi_new,old.eps_p_M, ...
    old.eps_p_AR);
[G_M0,G_AR0]=transfer_inventory(xi_old,xi_new,old.gamma_M, ...
    old.gamma_AR);
[D_M0,D_AR0]=transfer_inventory(xi_old,xi_new,old.rho_M,old.rho_AR);

xi_M_acc=old.xi_M_acc+abs(dxi);
mu=min(closure.chi_mu0*xi_new,ar_new);
free_zone=max(ar_new-mu,0);
if xi_new>tol
    chi_mu=mu/xi_new;
else
    chi_mu=closure.chi_mu0;
end
chi_f=closure.chi_f0_MPa*(exp(xi_new)-1);
chi_p=closure.chi_p_sat_MPa*(1-exp(-closure.b*xi_M_acc));
chi=chi_f+chi_p;
[mt_orientation,mt_orientation_owner]=transformation_orientation(closure);
beta=chi*mt_orientation;
chi_expected=closure.chi_f0_MPa*(exp(xi_new)-1) + ...
    closure.chi_p_sat_MPa*(1-exp(-closure.b*xi_M_acc));
chi_formula_residual=chi-chi_expected;
beta_formula_residual=beta-chi*mt_orientation;
assert(abs(chi_formula_residual)<=2e-12 && abs(beta_formula_residual)<=2e-12, ...
    'sma:candidate_arp_trip_inheritance_1d:BetaState', ...
    'Accepted beta state violates its declared chi/beta decomposition.');

if ar_new>tol
    if nargin<6 || isempty(guess), guess=[0;0]; end
    slip=double(guess(:));
    assert(numel(slip)==2 && all(isfinite(slip)), ...
        'sma:candidate_arp_trip_inheritance_1d:Guess', ...
        'Initial guess must be [Delta gamma_arp, Delta gamma_trip].');
    lower=[0;0]; upper=closure.max_local_slip_increment*[1;1];
    if ~closure.allow_arp, upper(1)=0; end
    if ~closure.allow_trip, upper(2)=0; end
    [rate_activity_arp,rate_activity_trip,rate_measure]= ...
        plastic_rate_activity(closure,free_zone,mu);
    if strcmp(rate_measure,'effective')
        % Conditional counterfactual: apply the power law to the
        % volume-weighted rate. Zhang Eqs. (40)--(41) instead use the
        % unbarred local rates; Eq. (14) weights their physical increments.
        % In this optional branch a vanishing active zone removes the
        % degenerate local-slip unknown. The default remains source-local.
        if rate_activity_arp<=tol, upper(1)=0; end
        if rate_activity_trip<=tol, upper(2)=0; end
    end
    slip=min(max(slip,lower),upper);
    active=find(upper>lower);
    fun=@(value)plastic_residual(value,G_AR0,ar_new,free_zone,mu, ...
        sigma_MPa,beta,dt,closure,lower,upper);
    converged=false; accepted_steps=0; rate_ratio_fallback_used=false;
    force_scalar_hardening=get_option(closure,'force_scalar_hardening',false);
    if force_scalar_hardening
        % Diagnostic/continuation path: solve the monotone shared-hardening
        % scalar root directly.  This is algebraically equivalent to the
        % two-slip projected system when both channels share tau_p, and is
        % useful for distinguishing a root-conditioning failure from a
        % constitutive incompatibility.  It is opt-in and does not alter
        % the historical default Newton path.
        [slip,q_override]=scalar_hardening_fallback(G_AR0,ar_new,free_zone,mu, ...
            sigma_MPa,beta,dt,closure,upper);
        [R,terms]=plastic_residual_from_q(slip,q_override,G_AR0,ar_new, ...
            free_zone,mu,sigma_MPa,beta,dt,closure,lower,upper);
        residual_norm=norm(R,inf);
        converged=residual_norm<=2e-10;
        rate_ratio_fallback_used=true;
        iteration=0;
    end
    for iteration=1:80
        if force_scalar_hardening, break; end
        R=fun(slip);
        residual_norm=norm(R,inf);
        if residual_norm<=2e-10
            converged=true;
            break;
        end
        J=sma.box_constrained_jacobian( ...
            @(value)active_output(fun,slip,active,value),slip(active), ...
            1e-7,lower(active),upper(active));
        R_active=R(active);
        if rcond(J)<1e-12
            step_active=-(J.'*J+1e-12*eye(numel(active)))\ ...
                (J.'*R_active);
        else
            step_active=-J\R_active;
        end
        accepted=false;
        for ls=1:24
            candidate=slip;
            candidate(active)=min(max(slip(active)+ ...
                (0.5^(ls-1))*step_active,lower(active)),upper(active));
            Rc=fun(candidate);
            if norm(Rc,inf)<residual_norm
                slip=candidate; accepted=true;
                accepted_steps=accepted_steps+1;
                break;
            end
        end
        if ~accepted, break; end
    end
    if ~force_scalar_hardening
        [R,terms]=fun(slip);
        residual_norm=norm(R,inf);
        converged=converged || residual_norm<=2e-10;
        if converged && plastic_rate_ratio_mismatch(slip,terms,dt,closure,lower,upper)
            % A tiny positive slip can satisfy the projected distance
            % residual while its direct 1/m power is an ill-conditioned
            % pseudo-root.  Re-evaluate this branch in finite q coordinates.
            converged=false;
        end
    end
    if ~converged
        % At activation corners mp=200 can make the two-flux Jacobian
        % effectively singular.  The common hardening stress is a monotone
        % scalar unknown, so use a bracketed log-domain fallback and still
        % accept only against the original projected root-form residual.
        [slip,q_override]=scalar_hardening_fallback(G_AR0,ar_new,free_zone,mu, ...
            sigma_MPa,beta,dt,closure,upper);
        [R,terms]=plastic_residual_from_q(slip,q_override,G_AR0,ar_new, ...
            free_zone,mu,sigma_MPa,beta,dt,closure,lower,upper);
        residual_norm=norm(R,inf);
        converged=residual_norm<=2e-10;
        rate_ratio_fallback_used=true;
    end
    assert(converged, ...
        'sma:candidate_arp_trip_inheritance_1d:Newton', ...
        ['ARP/TRIP shadow Newton failed; residual %.6g; slip=[%.6g %.6g]; ' ...
         'drive=[%.6g %.6g] MPa; rate_root=[%.6g %.6g] MPa; ' ...
         'pi=[%.6g %.6g] MPa.'], residual_norm, slip(1), slip(2), ...
        terms.drive_arp_MPa,terms.drive_trip_MPa,terms.rate_root_arp_MPa, ...
        terms.rate_root_trip_MPa,terms.pi_arp_MPa,terms.pi_trip_MPa);
else
    slip=[0;0]; R=[0;0]; residual_norm=0; iteration=0;
    accepted_steps=0; rate_ratio_fallback_used=false;
    terms=plastic_terms(slip,G_AR0,max(ar_new,eps),free_zone,mu, ...
        sigma_MPa,beta,dt,closure);
end

delta_gamma_arp=slip(1); delta_gamma_trip=slip(2);
N_arp=unit_direction(sigma_MPa,old.eps_p_AR);
N_trip=unit_direction(sigma_MPa+beta,mt_orientation);
delta_eps_arp=free_zone*delta_gamma_arp*N_arp;
delta_eps_trip=mu*delta_gamma_trip*N_trip;
P_AR=P_AR0+delta_eps_arp+delta_eps_trip;
G_AR=G_AR0+free_zone*delta_gamma_arp+mu*delta_gamma_trip;

% Kocks--Mecking law integrated analytically in q=sqrt(rho_AR), using
% the intrinsic average plastic increment in the current A+R inventory.
if ar_new>tol
    rho_AR_base=D_AR0/ar_new;
    intrinsic_average=(free_zone*delta_gamma_arp+ ...
        mu*delta_gamma_trip)/ar_new;
    q_old=sqrt(max(rho_AR_base,0));
    if closure.k2>0
        q_sat=closure.k1_per_m/closure.k2;
        q_new=q_sat+(q_old-q_sat)*exp(-0.5*closure.k2*intrinsic_average);
    else
        q_new=q_old+0.5*closure.k1_per_m*intrinsic_average;
    end
    rho_AR=max(q_new,0)^2;
    D_AR=ar_new*rho_AR;
    eps_p_AR=P_AR/ar_new;
    gamma_AR=G_AR/ar_new;
else
    rho_AR=old.rho_AR; D_AR=0;
    eps_p_AR=old.eps_p_AR; gamma_AR=old.gamma_AR;
    intrinsic_average=0;
end
if xi_new>tol
    rho_M=D_M0/xi_new; eps_p_M=P_M0/xi_new; gamma_M=G_M0/xi_new;
else
    rho_M=old.rho_M; eps_p_M=old.eps_p_M; gamma_M=old.gamma_M;
end

state=old;
state.xi_M=xi_new; state.xi_M_acc=xi_M_acc;
state.gamma_AR=max(gamma_AR,0); state.gamma_M=max(gamma_M,0);
state.eps_p_AR=eps_p_AR; state.eps_p_M=eps_p_M;
state.rho_AR=max(rho_AR,0); state.rho_M=max(rho_M,0);
state.eps_ARP=old.eps_ARP+delta_eps_arp;
state.eps_TRIP=old.eps_TRIP+delta_eps_trip;
state.chi=chi; state.beta=beta;

% Zhang Eqs. (53)--(60), driven by the unique conservative rho_AR/rho_M
% state above.  These constants are part of the source-owned plastic
% closure.  The defaults preserve the historical candidate exactly, while
% allowing the accepted-trial source card to bind the long-memory law to the
% same transaction instead of hiding it in this function.
a1=get_option(closure,'zhang_a1',0.6);
a2=get_option(closure,'zhang_a2',0.1);
cOmega=get_option(closure,'zhang_cOmega_MPa',92.0);
cxi=get_option(closure,'zhang_cxi_sqrt_m',3.5e-8);
etaM_N_per_K=get_option(closure,'zhang_etaM_N_per_K',1.23e-9);
etaA_N_per_K=get_option(closure,'zhang_etaA_N_per_K',1.0e-9);
cY=get_option(closure,'zhang_cY_MPa',4.5);
residual_event_policy=char(get_option(closure,'residual_m_event_policy', ...
    'clip_legacy'));
assert(any(strcmpi(residual_event_policy,{'clip_legacy','reject_crossing','event_boundary'})), ...
    'sma:candidate_arp_trip_inheritance_1d:ResidualEventPolicy', ...
    'Unknown residual-M event policy: %s.',residual_event_policy);
assert(all(isfinite([a1,a2,cOmega,cxi,etaM_N_per_K,etaA_N_per_K,cY])) && ...
    a1>=0 && a2>=0 && cOmega>=0 && cxi>=0 && ...
    etaM_N_per_K>=0 && etaA_N_per_K>=0 && cY>=0, ...
    'sma:candidate_arp_trip_inheritance_1d:HistoryConstants', ...
    'Zhang history constants must be finite nonnegative closure values.');
% The source eta coefficients produce Pa/K and are converted explicitly to
% MPa/K before storage and use in the MT free energy.
dLambda_M=abs(dxi);
xi_target_unclipped=cxi*sqrt(state.rho_M);
residual_event_projection_limit=0;
if any(strcmp(residual_event_policy,{'reject_crossing','event_boundary'}))
    % Eq.(53-b) uses the unbounded asymptotic target F.  Only the moving
    % admissibility boundary xi_Mr<=xi_M may reject a trial; clipping F to
    % xi_M changes the source ODE before the boundary event and is therefore
    % forbidden in this opt-in diagnostic mode.
    xi_target=xi_target_unclipped;
    xi_Mr_unclipped=xi_target+(old.xi_Mr-xi_target)*exp(-a1*dLambda_M);
    if strcmp(residual_event_policy,'reject_crossing')
        assert(xi_Mr_unclipped>=-tol && xi_Mr_unclipped<=xi_new+tol, ...
            'sma:candidate_arp_trip_inheritance_1d:ResidualEventCrossing', ...
            ['Unclipped Eq.(53-b) residual state crosses the moving boundary: ', ...
            'xi_Mr=%.17g, xi_M=%.17g.'],xi_Mr_unclipped,xi_new);
    else
        boundary_tol=get_option(closure,'residual_event_boundary_tolerance',2e-9);
        assert(isscalar(boundary_tol) && isfinite(boundary_tol) && boundary_tol>=0, ...
            'sma:candidate_arp_trip_inheritance_1d:ResidualEventTolerance', ...
            'Event-boundary tolerance must be finite and nonnegative.');
        residual_event_projection_limit=min(boundary_tol, ...
            256*eps(max([1,abs(xi_new),abs(xi_Mr_unclipped)])));
        assert(xi_Mr_unclipped>=-tol && ...
            xi_Mr_unclipped<=xi_new+residual_event_projection_limit, ...
            'sma:candidate_arp_trip_inheritance_1d:ResidualEventBoundaryMismatch', ...
            ['Event-boundary trial crosses xi_Mr=xi_M: ', ...
            'xi_Mr=%.17g, xi_M=%.17g, mismatch=%.17g.'], ...
            xi_Mr_unclipped,xi_new,xi_Mr_unclipped-xi_new);
    end
    state.xi_Mr=xi_Mr_unclipped;
else
    xi_target=min(max(xi_target_unclipped,0),xi_new);
    state.xi_Mr=xi_target+(old.xi_Mr-xi_target)*exp(-a1*dLambda_M);
    state.xi_Mr=min(max(state.xi_Mr,0),xi_new);
end
residual_event_boundary_projection_active=false;
residual_event_boundary_projection_amount=0;
if strcmp(residual_event_policy,'event_boundary')
    if state.xi_Mr>xi_new && ...
            state.xi_Mr<=xi_new+residual_event_projection_limit
        % Equation-level complementarity projection at the explicitly
        % derived xi_Mr=xi_M boundary.  This is not an experimental
        % endpoint assignment; the unclipped map and projection magnitude
        % remain recorded in the audit for restart/rollback verification.
        residual_event_boundary_projection_active=true;
        residual_event_boundary_projection_amount=state.xi_Mr-xi_new;
        state.xi_Mr=xi_new;
    end
end
Omega_expected=cOmega+(old.Omega-cOmega)*exp(-a2*dLambda_M);
state.Omega=Omega_expected;
state.Yp=cY+(old.Yp-cY)*exp(-a2*dLambda_M);
state.DeltaEtaP_MPaK=old.DeltaEtaP_MPaK+ ...
    (etaM_N_per_K*(state.rho_M-old.rho_M)- ...
    etaA_N_per_K*(state.rho_AR-old.rho_AR))/1e6;
state.Lambda_M=old.Lambda_M+dLambda_M;

plastic_inventory=xi_new*state.eps_p_M+ar_new*state.eps_p_AR;
slip_inventory=xi_new*state.gamma_M+ar_new*state.gamma_AR;
rho_inventory=xi_new*state.rho_M+ar_new*state.rho_AR;
old_plastic_inventory=xi_old*old.eps_p_M+ar_old*old.eps_p_AR;
old_slip_inventory=xi_old*old.gamma_M+ar_old*old.gamma_AR;
old_rho_inventory=xi_old*old.rho_M+ar_old*old.rho_AR;
rho_generation=D_AR-D_AR0;
audit=struct('delta_xi_M',dxi,'dxi_plus',dxi_plus, ...
    'dxi_minus',dxi_minus,'mu',mu,'free_zone',free_zone, ...
    'chi_mu',chi_mu,'chi_f_MPa',chi_f,'chi_p_MPa',chi_p, ...
    'chi_MPa',chi,'beta_MPa',beta,'N_arp',N_arp,'N_trip',N_trip, ...
    'beta_old_MPa',old.beta,'chi_old_MPa',old.chi, ...
    'chi_formula_residual_MPa',chi_formula_residual, ...
    'beta_formula_residual_MPa',beta_formula_residual, ...
    'beta_direction_MT',mt_orientation, ...
    'beta_state_convention','chi_f_plus_chi_p_then_beta_equals_chi_direction', ...
    'beta_intrinsic_path_channel_separate',true, ...
    'delta_gamma_arp_local',delta_gamma_arp, ...
    'delta_gamma_trip_local',delta_gamma_trip, ...
    'delta_eps_ARP',delta_eps_arp,'delta_eps_TRIP',delta_eps_trip, ...
    'intrinsic_average_slip',intrinsic_average, ...
    'plastic_inventory',plastic_inventory, ...
    'plastic_conservation_error',plastic_inventory- ...
        (old_plastic_inventory+delta_eps_arp+delta_eps_trip), ...
    'slip_inventory',slip_inventory, ...
    'slip_conservation_error',slip_inventory- ...
        (old_slip_inventory+free_zone*delta_gamma_arp+ ...
        mu*delta_gamma_trip), ...
    'rho_inventory',rho_inventory, ...
    'rho_balance_error',rho_inventory-(old_rho_inventory+rho_generation), ...
    'rho_generation',rho_generation,'residual',R, ...
    'xi_Mr',state.xi_Mr, ...
    'xi_Mr_unclipped',xi_target_unclipped+(old.xi_Mr-xi_target_unclipped)*exp(-a1*dLambda_M), ...
    'xi_Mr_target_unclipped',xi_target_unclipped, ...
    'xi_Mr_target_clipped',xi_target, ...
    'xi_Mr_target_clip_active',abs(xi_target-xi_target_unclipped)>100*eps, ...
    'xi_Mr_upper_clamp_active',abs(state.xi_Mr-xi_new)<=100*eps && ...
        xi_target_unclipped>=xi_new, ...
    'Omega_MPa',state.Omega, ...
    'Omega_old_MPa',old.Omega,'Omega_target_MPa',cOmega, ...
    'Omega_rate_per_abs_dM',a2, ...
    'Omega_delta_abs_M',dLambda_M, ...
    'Omega_evolution_increment_MPa',state.Omega-old.Omega, ...
    'Omega_evolution_residual_MPa',state.Omega-Omega_expected, ...
    'Omega_update_convention','exact_exponential_in_abs_delta_M', ...
    'Omega_formula_id','Zhang_Eq53a_Eq57_constant_target', ...
    'Omega_parameter_owner','candidate_arp_trip_inheritance_1d', ...
    'zhang_a1',a1,'zhang_a2',a2,'zhang_cOmega_MPa',cOmega, ...
    'zhang_cxi_sqrt_m',cxi,'zhang_etaM_N_per_K',etaM_N_per_K, ...
    'zhang_etaA_N_per_K',etaA_N_per_K,'zhang_cY_MPa',cY, ...
    'residual_m_event_policy',residual_event_policy, ...
    'residual_event_boundary_projection_active',residual_event_boundary_projection_active, ...
    'residual_event_boundary_projection_amount',residual_event_boundary_projection_amount, ...
    'residual_event_boundary_projection_limit',residual_event_projection_limit, ...
    'Yp_MPa',state.Yp,'DeltaEtaP_MPaK',state.DeltaEtaP_MPaK, ...
    'delta_Lambda_M',dLambda_M, ...
    'local_residual_norm',residual_norm,'iterations',iteration, ...
    'rate_ratio_fallback_used',rate_ratio_fallback_used, ...
    'accepted_steps',accepted_steps,'terms',terms, ...
    'observation_accessed',false,'production_connected',false);
if isfield(closure,'mt_transformation_orientation')
    audit.mt_transformation_orientation=mt_orientation;
    audit.mt_orientation_owner=mt_orientation_owner;
    audit.legacy_direction_MT_input=closure.direction_MT;
    audit.beta_state_convention= ...
        'chi_f_plus_chi_p_then_beta_equals_chi_transformation_orientation';
end
assert(abs(audit.plastic_conservation_error)<=2e-13 && ...
    abs(audit.slip_conservation_error)<=2e-13 && ...
    abs(audit.rho_balance_error)<=max(1,abs(rho_inventory))*2e-13, ...
    'sma:candidate_arp_trip_inheritance_1d:Conservation', ...
    'Accepted plastic/inheritance increment violates inventory balance.');
assert(state.xi_Mr>=-tol && state.xi_Mr<=state.xi_M+tol && ...
    state.Omega>=-tol && state.Yp>=-tol && state.Lambda_M>=old.Lambda_M-tol, ...
    'sma:candidate_arp_trip_inheritance_1d:HistoryBounds', ...
    'Accepted direct Zhang history violates its physical bounds.');
end

function [M_new,AR_new]=transfer_inventory(xi_old,xi_new,M_old,AR_old)
ar_old=1-xi_old; ar_new=1-xi_new;
M_inventory=xi_old*M_old; AR_inventory=ar_old*AR_old;
if xi_new>=xi_old
    transfer=(xi_new-xi_old)*AR_old;
    M_new=M_inventory+transfer;
    AR_new=AR_inventory-transfer;
else
    transfer=(xi_old-xi_new)*M_old;
    M_new=M_inventory-transfer;
    AR_new=AR_inventory+transfer;
end
% Remove roundoff only; no post-solve state correction is performed.
if abs(M_new)<100*eps, M_new=0; end
if abs(AR_new)<100*eps, AR_new=0; end
if ar_new==0, AR_new=0; end
end

function [R,terms]=plastic_residual(slip,G_AR0,ar_new,free_zone,mu, ...
    sigma,beta,dt,c,lower,upper)
terms=plastic_terms(slip,G_AR0,ar_new,free_zone,mu,sigma,beta,dt,c);
scale=max([abs(sigma),abs(beta),c.tau0_MPa,1]);
raw=[terms.pi_arp_MPa;terms.pi_trip_MPa]/scale;
R=slip-min(max(slip+raw,lower),upper);
end

function [R,terms]=plastic_residual_from_q(slip,q_override,G_AR0,ar_new, ...
    free_zone,mu,sigma,beta,dt,c,lower,upper)
% Preserve the finite rate-ratio root when mp=200 makes physical slip
% smaller than realmin.  The physical increment may remain zero in double
% precision while its force contribution is still represented by q.
terms=plastic_terms(slip,G_AR0,ar_new,free_zone,mu,sigma,beta,dt,c,q_override);
scale=max([abs(sigma),abs(beta),c.tau0_MPa,1]);
raw=[terms.pi_arp_MPa;terms.pi_trip_MPa]/scale;
R=slip-min(max(slip+raw,lower),upper);
end

function mismatch=plastic_rate_ratio_mismatch(slip,terms,dt,c,lower,upper)
% Detect a projected pseudo-root caused by the ill-conditioned direct
% DeltaGamma^(1/mp) evaluation.  Interior channels must have q_slip equal
% to the force-implied q_drive; true lower/upper KKT channels are exempt.
q_tol=get_option(c,'plastic_rate_ratio_tolerance',1e-8);
r0=get_option(c,'plastic_reference_rate_per_s',1);
drive=[terms.drive_arp_MPa,terms.drive_trip_MPa];
activity=[terms.rate_activity_arp,terms.rate_activity_trip];
mismatch=false;
for k=1:2
    if upper(k)<=lower(k), continue; end
    bound_tol=100*eps*max(1,upper(k));
    if upper(k)-slip(k)<=bound_tol, continue; end
    if slip(k)<=lower(k)+bound_tol
        q_slip=0;
    else
        q_slip=(activity(k)*slip(k)/(dt*r0))^(1/c.mp);
    end
    q_drive=max(drive(k)/c.tau0_MPa,0);
    if abs(q_slip-q_drive)>q_tol
        mismatch=true;
        return;
    end
end
end

function terms=plastic_terms(slip,G_AR0,ar_new,free_zone,mu, ...
    sigma,beta,dt,c,q_override)
gamma_AR=(G_AR0+free_zone*slip(1)+mu*slip(2))/ar_new;
tau_p=c.Hp_MPa*gamma_AR;
r0=get_option(c,'plastic_reference_rate_per_s',1);
[activity_arp,activity_trip,rate_measure]= ...
    plastic_rate_activity(c,free_zone,mu);
effective_increment_arp=activity_arp*slip(1);
effective_increment_trip=activity_trip*slip(2);
rate_ratio_override_used=nargin>=10 && ~isempty(q_override);
if rate_ratio_override_used
    q_override=q_override(:);
    assert(numel(q_override)==2 && all(isfinite(q_override)) && all(q_override>=0), ...
        'Rate-ratio override must contain two finite nonnegative values.');
    root_arp=c.tau0_MPa*q_override(1);
    root_trip=c.tau0_MPa*q_override(2);
else
    root_arp=c.tau0_MPa*(effective_increment_arp/(dt*r0))^(1/c.mp);
    root_trip=c.tau0_MPa*(effective_increment_trip/(dt*r0))^(1/c.mp);
end
% Keep the normal direction used by the signed plastic-strain update in
% the same transaction as the scalar driving force.  For nonzero stress
% this reproduces abs(sigma) algebraically; at zero stress the explicit
% direction policy provides the only auditable normal owner.
N_arp=unit_direction(sigma,1);
[mt_orientation,~]=transformation_orientation(c);
N_trip=unit_direction(sigma+beta,mt_orientation);
signed_drive_arp=sigma*N_arp-tau_p;
signed_drive_trip=(sigma+beta)*N_trip-tau_p;
drive_arp=abs(sigma)-tau_p;
drive_trip=abs(sigma+beta)-tau_p;
assert(abs(signed_drive_arp-drive_arp)<=1e-12*max(1,abs(drive_arp)) && ...
    abs(signed_drive_trip-drive_trip)<=1e-12*max(1,abs(drive_trip)), ...
    'Signed plastic normal and scalar drive are inconsistent.');
terms=struct('gamma_AR',gamma_AR,'tau_p_MPa',tau_p, ...
    'drive_arp_MPa',drive_arp,'drive_trip_MPa',drive_trip, ...
    'signed_drive_arp_MPa',signed_drive_arp, ...
    'signed_drive_trip_MPa',signed_drive_trip, ...
    'N_arp',N_arp,'N_trip',N_trip, ...
    'rate_root_arp_MPa',root_arp,'rate_root_trip_MPa',root_trip, ...
    'rate_measure',rate_measure,'rate_activity_arp',activity_arp, ...
    'rate_activity_trip',activity_trip, ...
    'rate_ratio_override_used',rate_ratio_override_used, ...
    'effective_increment_arp',effective_increment_arp, ...
    'effective_increment_trip',effective_increment_trip, ...
    'pi_arp_MPa',drive_arp-root_arp, ...
    'pi_trip_MPa',drive_trip-root_trip);
end

function output=active_output(fun,full_value,active,value)
full_value(active)=value;
output=fun(full_value);
output=output(active);
end

function [slip,q]=scalar_hardening_fallback(G_AR0,ar_new,free_zone,mu, ...
    sigma,beta,dt,c,upper)
[activity_arp,activity_trip]=plastic_rate_activity(c,free_zone,mu);
tau_base=c.Hp_MPa*G_AR0/ar_new;
tau_low=tau_base;
tau_high=max([abs(sigma),abs(sigma+beta),tau_base])+ ...
    c.Hp_MPa*(free_zone*upper(1)+mu*upper(2))/ar_new+c.tau0_MPa;
for k=1:200
    tau_mid=0.5*(tau_low+tau_high);
    [candidate1,~]=stable_slip_q(abs(sigma)-tau_mid,dt,c,upper(1),activity_arp);
    [candidate2,~]=stable_slip_q(abs(sigma+beta)-tau_mid,dt,c,upper(2),activity_trip);
    candidate=[candidate1;candidate2];
    balance=tau_mid-c.Hp_MPa*(G_AR0+free_zone*candidate(1)+ ...
        mu*candidate(2))/ar_new;
    if balance>0
        tau_high=tau_mid;
    else
        tau_low=tau_mid;
    end
end
tau=0.5*(tau_low+tau_high);
[slip1,q1]=stable_slip_q(abs(sigma)-tau,dt,c,upper(1),activity_arp);
[slip2,q2]=stable_slip_q(abs(sigma+beta)-tau,dt,c,upper(2),activity_trip);
slip=[slip1;slip2]; q=[q1;q2];
end

function [value,q]=stable_slip_q(overstress,dt,c,upper,activity)
if overstress<=0 || upper==0 || activity<=0
    value=0; q=0;
    return;
end
reference_rate=get_option(c,'plastic_reference_rate_per_s',1);
assert(isfinite(reference_rate) && reference_rate>0, ...
    'Plastic reference rate must be positive and finite.');
q=overstress/c.tau0_MPa;
log_value=log(dt*reference_rate/activity)+c.mp*log(q);
if log_value>=log(upper)
    value=upper;
    q=(activity*upper/(dt*reference_rate))^(1/c.mp);
elseif log_value<log(realmin)
    value=0;
else
    value=exp(log_value);
end
end

function [arp,trip,measure]=plastic_rate_activity(c,free_zone,mu)
measure=lower(char(get_option(c,'plastic_rate_measure','local')));
switch measure
    case 'local'
        arp=1; trip=1;
    case 'effective'
        arp=max(free_zone,0); trip=max(mu,0);
    otherwise
        error('sma:candidate_arp_trip_inheritance_1d:RateMeasure', ...
            'Unknown plastic rate measure: %s.',measure);
end
end

function value=get_option(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end

function [orientation,owner]=transformation_orientation(c)
% Eq.(25) orientation is distinct from the sign of the MT reaction. An
% explicit caller-owned orientation survives a protocol reversal. Omission
% preserves the legacy branch and its audit fields exactly.
if isfield(c,'mt_transformation_orientation')
    orientation=c.mt_transformation_orientation;
    assert(isnumeric(orientation) && isreal(orientation) && ...
        isscalar(orientation) && isfinite(orientation) && any(orientation==[-1,1]), ...
        'sma:candidate_arp_trip_inheritance_1d:MTOrientation', ...
        'Explicit MT transformation orientation must be scalar +1 or -1.');
    assert(isfield(c,'mt_orientation_owner') && ...
        (ischar(c.mt_orientation_owner) || isstring(c.mt_orientation_owner)), ...
        'sma:candidate_arp_trip_inheritance_1d:MTOrientationOwner', ...
        'Explicit MT orientation requires its caller-owned provenance.');
    owner=char(c.mt_orientation_owner);
    assert(isrow(owner) && ~isempty(strtrim(owner)), ...
        'sma:candidate_arp_trip_inheritance_1d:MTOrientationOwner', ...
        'MT orientation owner must be one nonempty text value.');
else
    assert(~isfield(c,'mt_orientation_owner'), ...
        'sma:candidate_arp_trip_inheritance_1d:MTOrientationOwner', ...
        'MT orientation provenance cannot be supplied without an orientation.');
    orientation=c.direction_MT;
    owner='legacy_direction_MT';
end
end

function direction=unit_direction(value,fallback)
direction=sign(value);
if direction==0, direction=sign(fallback); end
if direction==0, direction=1; end
end
