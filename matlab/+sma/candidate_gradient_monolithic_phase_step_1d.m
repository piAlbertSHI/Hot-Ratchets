function [phase,audit]=candidate_gradient_monolithic_phase_step_1d( ...
    p,old,macro_strain,T,dt,closure,guess)
%SMA.CANDIDATE_GRADIENT_MONOLITHIC_PHASE_STEP_1D Coupled spatial RT/MT solve.
% Unknowns are one common stress and every local RT/MT increment. This is
% the phase-only inner solve required before the ARP/TRIP consistency loop.

if nargin<7 || isempty(guess), guess=[]; end
n=numel(old); assert(iscell(old) && n>=1,'At least one phase cell is required.');
T=expand(T,n); w=expand(closure.weights,n); resistance=expand(closure.resistance_shift_MPa,n);
if ~isfield(closure,'enable_path_accommodation')
    closure.enable_path_accommodation=false;
end
assert(isscalar(closure.enable_path_accommodation) && ...
    any(closure.enable_path_accommodation==[false,true]), ...
    'The path-accommodation switch must be a scalar logical value.');
if closure.enable_path_accommodation
    required_path={'path_accommodation_moduli_MPa', ...
        'path_accommodation_resistances_MPa'};
    for k=1:numel(required_path)
        assert(isfield(closure,required_path{k}), ...
            'Missing path-accommodation closure field %s.',required_path{k});
    end
    nh=numel(closure.path_accommodation_moduli_MPa);
    assert(nh>=1 && ...
        numel(closure.path_accommodation_resistances_MPa)==nh, ...
        'Path-accommodation modulus and resistance widths must match.');
    for i=1:n
        assert(isfield(old{i},'path_accommodation_strain') && ...
            numel(old{i}.path_accommodation_strain)==nh, ...
            'Every committed phase cell must carry the path-accommodation state.');
    end
end
if ~isfield(closure,'enable_thermal_accommodation')
    closure.enable_thermal_accommodation=false;
end
assert(isscalar(closure.enable_thermal_accommodation) && ...
    any(closure.enable_thermal_accommodation==[false,true]), ...
    'The thermal-accommodation switch must be a scalar logical value.');
if closure.enable_thermal_accommodation
    required_thermal={'thermal_accommodation_moduli_MPa', ...
        'thermal_accommodation_resistances_MPa','committed_temperature_K'};
    for k=1:numel(required_thermal)
        assert(isfield(closure,required_thermal{k}), ...
            'Missing thermal-accommodation closure field %s.',required_thermal{k});
    end
    nt=numel(closure.thermal_accommodation_moduli_MPa);
    assert(nt>=1 && ...
        numel(closure.thermal_accommodation_resistances_MPa)==nt, ...
        'Thermal-accommodation modulus and resistance widths must match.');
    closure.committed_temperature_K=expand(closure.committed_temperature_K,n);
    for i=1:n
        assert(isfield(old{i},'thermal_accommodation_strain') && ...
            numel(old{i}.thermal_accommodation_strain)==nt, ...
            'Every committed phase cell must carry the thermal-accommodation state.');
    end
end
% Zhang's published kinetic law is written in terms of a dimensionless
% rate ratio, but the implementation must still choose the time scale used
% to convert a finite increment into a rate.  Keep the historical values as
% the strict default, while allowing an explicit diagnostic continuation to
% expose this otherwise implicit normalization.
kinetic_exponent=get_option(closure,'transformation_exponent',200);
kinetic_exponent_R=get_option(closure,'transformation_exponent_R',kinetic_exponent);
kinetic_exponent_M=get_option(closure,'transformation_exponent_M',kinetic_exponent);
kinetic_reference_rate=get_option(closure,'transformation_reference_rate_per_s',1);
assert(isfinite(kinetic_exponent_R) && kinetic_exponent_R>0 && ...
    isfinite(kinetic_exponent_M) && kinetic_exponent_M>0 && ...
    isfinite(kinetic_reference_rate) && kinetic_reference_rate>0, ...
    'Transformation kinetic normalization must be positive and finite.');
assert(all(w>0) && abs(sum(w)-1)<=1e-12 && abs(w.'*resistance)<=1e-12*max(1,max(abs(resistance))), ...
    'Positive normalized weights and zero-mean resistance are required.');
population_mode=isfield(closure,'population_count');
if population_mode && closure.kappa_MPa_length2~=0
    np=round(closure.population_count); assert(np>=1 && mod(n,np)==0, ...
        'Population count must divide the flattened state width.');
    nx=n/np; pop_weights=expand(closure.population_weights,np);
    [Dxx,~,operator]=sma.candidate_neumann_laplacian_1d(nx,closure.spatial_dx);
elseif closure.kappa_MPa_length2==0
    % Exact neutral branch: avoid introducing a different sparse-matrix
    % evaluation at ell=0, where population-space and legacy kernels must
    % be algebraically identical.
    if population_mode
        np=max(1,round(closure.population_count)); nx=n/np;
        pop_weights=expand(closure.population_weights,np);
    else
        np=1; nx=n; pop_weights=1;
    end
    Dxx=[]; operator=struct('type','zero_gradient');
else
    np=1; nx=n; pop_weights=1;
    [Dxx,~,operator]=sma.candidate_neumann_laplacian_1d(n,closure.dx);
end
xi0=cellfun(@(s)s.xi_M,old(:));
xi_res=expand(closure.xi_res,n); xi_recent=expand(closure.xi_recent,n);
Btr=expand(closure.Btr_MPa,n); deta_p=expand(closure.delta_eta_p_MPaK,n);
Ybase=expand(closure.Y_M_base_MPa,n); capacity=expand(closure.f_transformable,n);
Ymult=expand(get_option(closure,'resistance_multiplier',1),n);
% This function later owns numeric arrays named lower/upper for the bounded
% solve.  Avoid calling the text function LOWER here because MATLAB resolves
% that identifier as the local array throughout the function scope.
source_force_representation=char(get_option(closure, ...
    'source_force_representation','main_eq36b'));
assert(any(strcmp(source_force_representation, ...
    {'main_eq36b','supplement_s4a'})), ...
    'Unknown Zhang source-force representation: %s.',source_force_representation);
phase_composition_mode=char(get_option(closure, ...
    'phase_composition_mode','origin_memory'));
if strcmpi(phase_composition_mode,'origin_memory')
    phase_composition_mode='origin_memory';
elseif strcmpi(phase_composition_mode,'source_current_mixture')
    phase_composition_mode='source_current_mixture';
else
    error('Unknown Zhang phase-composition mode: %s.',phase_composition_mode);
end
closure.phase_composition_mode=phase_composition_mode;
Omega_state=zeros(n,1);
source_receipt_required=get_option(closure,'require_zhang_source_receipt',false);
source_branch_envelope=get_option(closure,'source_branch_envelope',struct());
branch_bound=isstruct(source_branch_envelope) && ...
    isfield(source_branch_envelope,'semantic_card_id');
if source_receipt_required
    eq39_ablation=get_option(closure,'eq39_ablation',canonical_eq39_ablation());
    if eq39_ablation.enabled
        assert(eq39_ablation.m_R==p.zhang_m_tr && ...
            eq39_ablation.reference_rate_per_s==p.zhang_reference_rate_per_s && ...
            ~eq39_ablation.source_authoritative && ...
            ~eq39_ablation.production_connected && ...
            ~eq39_ablation.experimental_acceptance && ...
            ~eq39_ablation.observation_accessed && ...
            ~eq39_ablation.experimental_backfill, ...
            'sma:candidate_gradient_monolithic_phase_step_1d:SourceReceipt', ...
            'Eq39 alternative-m branch is validation-only and must remain fail-closed.');
    else
        eq39_ablation=canonical_eq39_ablation();
    end
    assert(isfield(closure,'Y_M_threshold_receipt') && ...
        iscell(closure.Y_M_threshold_receipt) && ...
        numel(closure.Y_M_threshold_receipt)==n, ...
        'Direct Zhang gradient mode requires one threshold receipt per cell.');
    for i=1:n
        receipt=closure.Y_M_threshold_receipt{i};
        required_receipt={'material_case_id','xi_recent','Yp_MPa','Yf_MPa', ...
            'Y_M_base_MPa','formula_id','source_role','source_card_id', ...
            'Yf0_MPa','d_MPa','n_xi','source_transformation_exponent', ...
            'source_reference_rate_per_s','exponent','reference_rate_per_s', ...
            'source'};
        if branch_bound
            required_receipt=[required_receipt, ...
                {'source_branch_id','source_semantic_card_id', ...
                 'source_branch_card_digest','source_branch_authority'}];
        end
        for k=1:numel(required_receipt)
            assert(isfield(receipt,required_receipt{k}), ...
                'sma:candidate_gradient_monolithic_phase_step_1d:SourceReceipt', ...
                'Gradient Zhang receipt is missing %s.',required_receipt{k});
        end
        if eq39_ablation.enabled
            expected_role='non_source_eq39_ablation';
            expected_source='validation_only_counterfactual_m_branch';
            expected_authoritative=false;
            expected_exponent=eq39_ablation.m_M;
        else
            expected_role='kinetic_threshold';
            expected_source='Zhang_Table2_Eqs54_56';
            expected_authoritative=true;
            expected_exponent=p.zhang_m_tr;
        end
        assert(receipt.material_case_id==1 && ...
            abs(receipt.xi_recent-xi_recent(i))<=2e-14 && ...
            receipt.Yp_MPa>=0 && receipt.Yf_MPa>0 && ...
            abs(receipt.Y_M_base_MPa-Ybase(i))<=2e-12 && ...
            abs(receipt.Y_M_base_MPa-(receipt.Yf_MPa-receipt.Yp_MPa))<=2e-12 && ...
            abs(receipt.Yf0_MPa-p.zhang_Yf0_MPa)<=2e-14 && ...
            abs(receipt.d_MPa-p.zhang_Y_d_MPa)<=2e-14 && ...
            abs(receipt.n_xi-p.zhang_Y_n_xi)<=2e-14 && ...
            abs(receipt.source_transformation_exponent-p.zhang_m_tr)<=2e-14 && ...
            abs(receipt.source_reference_rate_per_s- ...
                p.zhang_reference_rate_per_s)<=2e-14 && ...
            strcmp(receipt.formula_id,'Zhang_Eq39_Yf_minus_Yp') && ...
            strcmp(receipt.source_role,expected_role) && ...
            strcmp(receipt.source_card_id,'niti_zhang2022_table2') && ...
            receipt.exponent==expected_exponent && ...
            receipt.reference_rate_per_s==p.zhang_reference_rate_per_s && ...
            receipt.exponent==kinetic_exponent_M && ...
            (eq39_ablation.enabled || receipt.exponent==kinetic_exponent_R) && ...
            receipt.reference_rate_per_s==kinetic_reference_rate && ...
            strcmp(receipt.source,expected_source) && ...
            isfield(receipt,'source_authoritative') && ...
            receipt.source_authoritative==expected_authoritative && ...
            isfield(receipt,'ablation_id') && ...
            strcmp(receipt.ablation_id,eq39_ablation.ablation_id) && ...
            isfield(receipt,'parent_source_receipt_digest') && ...
            strcmp(receipt.parent_source_receipt_digest,eq39_ablation.parent_source_receipt_digest), ...
            'sma:candidate_gradient_monolithic_phase_step_1d:SourceReceipt', ...
            ['Gradient Zhang threshold/kinetic receipt is inconsistent ' ...
             'with the executed closure.']);
        if branch_bound
            assert(strcmp(receipt.source_branch_id,source_branch_envelope.branch_id) && ...
                strcmp(receipt.source_semantic_card_id,source_branch_envelope.semantic_card_id) && ...
                strcmp(receipt.source_branch_card_digest,source_branch_envelope.card_payload_digest) && ...
                strcmp(receipt.source_branch_authority,source_branch_envelope.authority_class), ...
                'sma:candidate_gradient_monolithic_phase_step_1d:SourceReceipt', ...
                'Gradient Zhang receipt is bound to the wrong semantic source branch.');
        end
    end
end
% The full gradient step must pass one explicit history packet into this
% kernel.  This is deliberately a wiring/identity check only: it does not
% manufacture any of the still-unresolved configurational, RT-common,
% Omega, or beta-companion force providers.  A missing or stale packet is
% rejected instead of allowing a formally converged residual to mix state
% values from different trials.
history_receipt_required=get_option(closure,'require_history_state_receipt',false);
if history_receipt_required
    assert(isfield(closure,'history_state_receipt') && ...
        iscell(closure.history_state_receipt) && ...
        numel(closure.history_state_receipt)==n, ...
        'Gradient history closure requires one state receipt per cell.');
    for i=1:n
        receipt=closure.history_state_receipt{i};
        required_history={'material_case_id','history_owner','history_state_token', ...
            'xi_Mr','Omega_MPa','Yp_MPa','DeltaEtaP_MPaK','Btr_MPa', ...
            'Yf_MPa','Y_M_base_MPa','source_card_id','formula_id'};
        if branch_bound
            required_history=[required_history, ...
                {'source_branch_id','source_semantic_card_id', ...
                 'source_branch_card_digest','source_branch_authority'}];
        end
        for k=1:numel(required_history)
            assert(isfield(receipt,required_history{k}), ...
                'sma:candidate_gradient_monolithic_phase_step_1d:HistoryReceipt', ...
                'Gradient history receipt is missing %s.',required_history{k});
        end
        assert(receipt.material_case_id==p.case_id && ...
            strcmp(receipt.history_owner,'candidate_gradient_monolithic_full_step_1d') && ...
            ~isempty(receipt.history_state_token) && ...
            isfinite(receipt.xi_Mr) && abs(receipt.xi_Mr-xi_res(i))<=2e-14 && ...
            isfinite(receipt.Omega_MPa) && isfinite(receipt.Yp_MPa) && ...
            isfinite(receipt.DeltaEtaP_MPaK) && ...
            abs(receipt.Btr_MPa-Btr(i))<=2e-12 && ...
            abs(receipt.Yf_MPa-receipt.Yp_MPa-Ybase(i))<=2e-12 && ...
            abs(receipt.Y_M_base_MPa-Ybase(i))<=2e-12 && ...
            strcmp(receipt.source_card_id,'niti_zhang2022_table2') && ...
            strcmp(receipt.formula_id,'Zhang_Eq39_Yf_minus_Yp') && ...
            strcmp(receipt.history_state_token,history_state_token(receipt)), ...
            'sma:candidate_gradient_monolithic_phase_step_1d:HistoryReceipt', ...
             'Gradient history receipt is inconsistent with the phase closure.');
        if branch_bound
            assert(strcmp(receipt.source_branch_id,source_branch_envelope.branch_id) && ...
                strcmp(receipt.source_semantic_card_id,source_branch_envelope.semantic_card_id) && ...
                strcmp(receipt.source_branch_card_digest,source_branch_envelope.card_payload_digest) && ...
                strcmp(receipt.source_branch_authority,source_branch_envelope.authority_class), ...
                'sma:candidate_gradient_monolithic_phase_step_1d:HistoryReceipt', ...
                'Gradient history receipt is bound to the wrong semantic source branch.');
        end
        Omega_state(i)=receipt.Omega_MPa;
    end
end
if strcmp(source_force_representation,'supplement_s4a')
    assert(history_receipt_required, ...
        'Supplementary S-4(a) force requires an exact committed Omega receipt.');
end
direction_R=expand(closure.direction_R,n); direction_M=expand(closure.direction_M,n);
% In committed-split mode xi_res is a committed internal fraction and must
% not exceed the committed xi_M.  A same-trial history fixed point instead
% supplies a trial iterate; its admissible upper bound is the transformable
% capacity, while the accepted plastic update still enforces xi_Mr<=xi_M.
coupling_mode=char(get_option(closure,'residual_martensite_coupling', ...
    'committed_split'));
if strcmpi(coupling_mode,'committed_split')
    xi_res_upper=xi0;
    residual_bound_tolerance=0;
else
    xi_res_upper=min(max(capacity,0),1);
    residual_bound_tolerance=100*eps;
end
valid_direction=all(direction_R==1 | direction_R==-1) && ...
    all(direction_M==1 | direction_M==-1);
valid_residual=all(xi_res>=0 & xi_res<=xi_res_upper+residual_bound_tolerance);
valid_capacity=all(capacity>=xi0 & capacity<=1);
valid_resistance=all(Ybase.*Ymult+resistance>0);
assert(valid_direction && valid_residual && valid_capacity && valid_resistance, ...
    ['Monolithic phase closure is inadmissible: direction=%d, ' ...
     'residual=%d, capacity=%d, resistance=%d, coupling=%s.'], ...
    valid_direction,valid_residual,valid_capacity,valid_resistance,coupling_mode);

% Validation-only residual-M event mode may supply an analytically derived
% lower bound for reverse M.  The bound is an admissibility surface from
% Eq.(53-b), not a post-solve clipping rule; the constitutive residual is
% still solved at the bound and the plastic kernel must close xi_Mr=xi_M.
has_event_lower=isfield(closure,'residual_event_xi_M_lower');
if has_event_lower
    event_xi_lower=expand(closure.residual_event_xi_M_lower,n);
    reverse_event_applicable=direction_M<0;
    % The event surface bounds reverse depletion only. In a forward
    % same-trial history solve Mr may grow above the old reverse-event
    % metadata, which is never used to constrain a positive M increment.
    % Retain finite/physical input checks and the original reverse bounds.
    assert(all(isfinite(event_xi_lower)) && ...
        all(event_xi_lower>=-2e-10 & event_xi_lower<=1+2e-10) && ...
        all(event_xi_lower(reverse_event_applicable)>=xi_res(reverse_event_applicable)-2e-10) && ...
        all(event_xi_lower(reverse_event_applicable)<=xi0(reverse_event_applicable)+2e-10), ...
        'sma:candidate_gradient_monolithic_phase_step_1d:ResidualEventLower', ...
        'Residual event lower bound is outside the committed phase interval.');
else
    event_xi_lower=xi_res;
end

lower=zeros(1+2*n,1); upper=zeros(1+2*n,1);
lower(1)=-p.numerics.stress_solve_bound; upper(1)=p.numerics.stress_solve_bound;
for i=1:n
    lower(1+i)=-old{i}.xi_R; upper(1+i)=old{i}.xi_A;
    lower(1+n+i)=-old{i}.xi_M+xi_res(i);
    if has_event_lower && direction_M(i)<0
        lower(1+n+i)=max(lower(1+n+i),event_xi_lower(i)-old{i}.xi_M);
    end
    upper(1+n+i)=min(old{i}.xi_A+old{i}.xi_R,capacity(i)-old{i}.xi_M);
    if direction_R(i)>0, lower(1+i)=0; else, upper(1+i)=0; end
    if direction_M(i)>0, lower(1+n+i)=0; else, upper(1+n+i)=0; end
end
if isempty(guess)
    y=zeros(1+2*n,1); y(1)=min(max(p.EA*macro_strain,lower(1)),upper(1));
else
    y=guess(:); assert(numel(y)==1+2*n && all(isfinite(y)), ...
        'The monolithic guess has an invalid width.');
end
% Opt-in force control replaces only the macro compatibility equation.
% The known stress is an external boundary value, never an observed state.
control_mode=char(get_option(closure,'candidate_mechanical_control','prescribed_strain'));
assert(any(strcmp(control_mode,{'prescribed_strain','validation_prescribed_stress'})), ...
    'Unknown candidate mechanical control owner.');
if strcmp(control_mode,'validation_prescribed_stress')
    assert(isfield(closure,'candidate_target_stress_MPa') && ...
        isscalar(closure.candidate_target_stress_MPa) && ...
        isfinite(closure.candidate_target_stress_MPa) && ...
        closure.candidate_target_stress_MPa>=lower(1) && ...
        closure.candidate_target_stress_MPa<=upper(1), ...
        'Explicit admissible prescribed stress is required.');
    y(1)=closure.candidate_target_stress_MPa;
end
y=min(max(y,lower),upper);
fun=@(v)evaluate(v,p,old,macro_strain,T,dt,closure,Dxx,xi0, ...
    xi_res,xi_recent,Btr,deta_p,Ybase,Ymult,resistance,capacity,direction_R,direction_M,w, ...
    population_mode,np,nx,pop_weights,kinetic_exponent_R,kinetic_exponent_M,kinetic_reference_rate, ...
    source_force_representation,Omega_state);
eval_args={p,old,macro_strain,T,dt,closure,Dxx,xi0,xi_res,xi_recent, ...
    Btr,deta_p,Ybase,Ymult,resistance,capacity,direction_R,direction_M,w, ...
    population_mode,np,nx,pop_weights,kinetic_exponent_R,kinetic_exponent_M, ...
    kinetic_reference_rate,source_force_representation,Omega_state};
% The projected complementarity residual is only piecewise smooth at
% activation/deactivation bounds. Central differences remain the historical
% default for contract compatibility; forward differences are available as
% an explicit candidate option when an activation kink is encountered.
solver_mode=get_option(closure,'phase_solver_mode','lsqnonlin');
transformed=false;
transform_kind='';
q_trial_final=[];
if strcmpi(solver_mode,'newton_q')
    % Solve in the natural rate-ratio coordinate q=(|dLambda|/(dt*r0))^(1/m).
    % In dLambda coordinates the derivative of the power law is singular at
    % activation (the low-amplitude root can be 1e-9 or smaller).  The q
    % coordinate makes the kinetic residual locally affine while retaining
    % exactly the same physical increment and admissible bounds.
    scale_dt=dt*kinetic_reference_rate;
    sgn=[direction_R(:);direction_M(:)];
    mcoord=[kinetic_exponent_R*ones(n,1);kinetic_exponent_M*ones(n,1)];
    lower_phys=lower; upper_phys=upper;
    seed_residual=fun(y);
    seed_scale=max(p.reference_stress,1);
    seed_resistance=[closure.Y_R_MPa*ones(n,1); ...
        Ybase.*Ymult+resistance];
    q_seed=abs(seed_residual(2:end))*seed_scale./seed_resistance;
    z=map_y_to_q(y,sgn,scale_dt,mcoord);
    lower=zeros(size(lower_phys)); upper=zeros(size(upper_phys));
    lower(1)=lower_phys(1); upper(1)=upper_phys(1);
    for kk=2:numel(y)
        lim=max(abs([lower_phys(kk),upper_phys(kk)])/scale_dt,0);
        lower(kk)=0; upper(kk)=max(lim)^(1/mcoord(kk-1));
        if z(kk)==0 && q_seed(kk-1)>0
            % The m>>1 rate law has a vanishing physical-increment
            % derivative at q=0. Seed q from the zero-increment driving
            % residual (q=|pi|/Y), then clamp to the admissible box. This
            % changes only the nonlinear initial guess, never the equation.
            z(kk)=min(q_seed(kk-1),upper(kk));
        end
    end
    q_seed_policy=char(get_option(closure, ...
        'phase_q_seed_policy','residual_scaled'));
    if strcmpi(q_seed_policy,'active_bound')
        % At a projected phase-capacity boundary, the zero-increment
        % residual is a phase-fraction distance, not a kinetic force ratio.
        % Scaling that distance by Y seeds q near zero; for m=200 the map
        % dLambda=dt*q^m then underflows and its Jacobian cannot reach the
        % active bound. Seed only nonzero projected components at their
        % admissible q bound. This changes the nonlinear initial iterate,
        % never the physical box, residual, or acceptance tolerance.
        seed_active_tolerance=get_option(closure, ...
            'phase_q_seed_active_tolerance',100*eps);
        for kk=2:numel(y)
            if abs(seed_residual(kk))>seed_active_tolerance
                z(kk)=upper(kk);
            end
        end
    elseif ~strcmpi(q_seed_policy,'residual_scaled')
        error('Unknown phase q seed policy %s.',q_seed_policy);
    end
    y=min(max(z,lower),upper); fun_phys=fun;
    if get_option(closure,'strict_source_kinetics',false)
        fun=@(v)evaluate_q_coordinate(v,'q',sgn,scale_dt,mcoord, ...
            lower_phys,upper_phys,eval_args);
    else
        fun=@(v)fun_phys(min(max(map_q_to_y(v,sgn,scale_dt,mcoord), ...
            lower_phys),upper_phys));
    end
    transformed=true; transform_kind='q';
elseif strcmpi(solver_mode,'newton_logq')
    % Logarithmic rate-ratio coordinate for very large m.  The q coordinate
    % compresses a many-decade physical increment range into a narrow band
    % near q=1 when m=200.  With u=log(q), the same constitutive law is
    % dLambda=sign*dt*r0*exp(m*u), while the kinetic resistance is linear in
    % exp(u).  The finite lower bound maps to realmin-scale increments and is
    % therefore the numerical representation of the exact inactive bound.
    scale_dt=dt*kinetic_reference_rate;
    sgn=[direction_R(:);direction_M(:)];
    mcoord=[kinetic_exponent_R*ones(n,1);kinetic_exponent_M*ones(n,1)];
    lower_phys=lower; upper_phys=upper;
    seed_residual=fun(y);
    seed_scale=max(p.reference_stress,1);
    seed_resistance=[closure.Y_R_MPa*ones(n,1); ...
        Ybase.*Ymult+resistance];
    q_seed=abs(seed_residual(2:end))*seed_scale./seed_resistance;
    z=map_y_to_logq(y,sgn,scale_dt,mcoord);
    lower=zeros(size(lower_phys)); upper=zeros(size(upper_phys));
    lower(1)=lower_phys(1); upper(1)=upper_phys(1);
    for kk=2:numel(y)
        lim=max(abs([lower_phys(kk),upper_phys(kk)])/scale_dt,0);
        upper(kk)=log(max(max(lim),realmin('double')))/mcoord(kk-1);
        lower(kk)=log(realmin('double'))/mcoord(kk-1);
        if q_seed(kk-1)>0
            z(kk)=max(z(kk),log(max(q_seed(kk-1),realmin('double'))));
        end
    end
    y=min(max(z,lower),upper);
    fun_phys=fun;
    if get_option(closure,'strict_source_kinetics',false)
        fun=@(v)evaluate_q_coordinate(v,'logq',sgn,scale_dt,mcoord, ...
            lower_phys,upper_phys,eval_args);
    else
        fun=@(v)fun_phys(min(max(map_logq_to_y(v,sgn,scale_dt,mcoord), ...
            lower_phys),upper_phys));
    end
    transformed=true; transform_kind='logq';
end
if strcmpi(solver_mode,'newton') || transformed
    % Toolbox-independent diagnostic semismooth Newton continuation.  The
    % bounded least-squares path remains the production default.  Every
    % trial is projected into the constitutive box and accepted only when
    % the original projected residual decreases; no endpoint state is set.
    maxit=get_option(closure,'phase_max_iterations',1000);
    maxfev=get_option(closure,'phase_max_function_evaluations',50000);
    residual=fun(y); fn=1;
    exitflag=0; output=struct('iterations',0,'funcCount',fn);
    for it=1:maxit
        output.iterations=it;
        if norm(residual,inf)<=get_option(closure,'residual_tolerance',2e-8)
            exitflag=1; break;
        end
        nv=numel(y); J=zeros(numel(residual),nv);
        for jj=1:nv
            % Resolve the very small finite increments generated by the
            % power-law rate near an activation surface.  A fixed 1e-7
            % increment is larger than the physical root at low amplitude
            % (and makes the projected Jacobian report a false plateau).
            if jj==1
                h=1e-7*max(1,abs(y(jj)));
            elseif transformed
                % q and log(q) are O(1) kinetic coordinates even when the
                % physical phase increment is zero.  The former formula
                % used max(|q|,1e-6), which produced h=1e-12 at activation
                % and numerically erased the rate-residual derivative.
                % Perturb the transformed coordinate itself at its natural
                % scale; the exact map below still evaluates the unchanged
                % physical constitutive equations.
                h=1e-6*max(1,abs(y(jj)));
            else
                h=1e-6*max(abs(y(jj)),1e-6);
            end
            yp=y; ym=y; yp(jj)=min(y(jj)+h,upper(jj)); ym(jj)=max(y(jj)-h,lower(jj));
            if yp(jj)==ym(jj), J(:,jj)=0; else
                rp=fun(yp); rm=fun(ym); fn=fn+2;
                J(:,jj)=(rp-rm)/(yp(jj)-ym(jj));
            end
            if fn>=maxfev, break; end
        end
        if any(~isfinite(J(:))), break; end
        % Use an SVD-based reciprocal condition estimate instead of
        % rcond(J): the latter emits a warning at every near-active-set
        % trial and can dominate a long thermal run's wall time/output.
        [U,S,V]=svd(J,'econ'); sj=diag(S);
        if isempty(sj) || sj(1)==0
            step=zeros(nv,1);
        elseif sj(end)/sj(1)<1e-12
            % Damped SVD solve (Tikhonov regularization) avoids the
            % near-singular backslash warning without changing the
            % constitutive residual or admissible box.
            lam=1e-12*sj(1)^2;
            step=-V*((sj./(sj.^2+lam)).*(U.'*residual));
        else
            step=-V*((U.'*residual)./sj);
        end
        accepted=false; oldnorm=norm(residual,inf);
        for ls=0:24
            trial=min(max(y+(0.5^ls)*step,lower),upper);
            rt=fun(trial); fn=fn+1;
            if norm(rt,inf)<oldnorm
                y=trial; residual=rt; accepted=true; break;
            end
            if fn>=maxfev, break; end
        end
        % At an activation bound the projected residual has a kink.  A
        % finite-difference Newton step can remain a few ulps inside the
        % admissible interval even though the complementarity root is the
        % bound itself.  Evaluate an explicit active-set projection as a
        % mathematically admissible trial (it is not an imposed endpoint
        % state): only variables already within the Jacobian stencil of a
        % bound are projected, and the trial is accepted only if the
        % original residual decreases.
        if ~accepted && fn<maxfev
            active_trial=y; hact=1e-7*max(1,abs(y));
            for kk=2:numel(y)
                if y(kk)-lower(kk)<=hact(kk), active_trial(kk)=lower(kk); end
                if upper(kk)-y(kk)<=hact(kk), active_trial(kk)=upper(kk); end
            end
            if any(active_trial~=y)
                ra=fun(active_trial); fn=fn+1;
                if norm(ra,inf)<oldnorm
                    y=active_trial; residual=ra; accepted=true;
                end
            end
        end
        output.funcCount=fn;
        if ~accepted || fn>=maxfev, break; end
    end
    output.funcCount=fn;
    % A semismooth Newton step can terminate at a projected active-set kink
    % with a residual only marginally above the strict contract tolerance.
    % Before declaring the constitutive transaction failed, perform one
    % bounded least-squares polish in the same coordinates and box. This is
    % an equation solve (no endpoint/observation data are available), and the
    % polished point is accepted only after re-evaluating the original
    % projected residual against the unchanged tolerance.
    strict_tol=get_option(closure,'phase_residual_tolerance',2e-8);
    if norm(residual,inf)>strict_tol && fn<maxfev
        try
            polish_opts=optimoptions('lsqnonlin','Display','off', ...
                'FunctionTolerance',1e-16,'StepTolerance',1e-15, ...
                'OptimalityTolerance',1e-14,'MaxIterations',500, ...
                'MaxFunctionEvaluations',max(2000,min(20000,maxfev-fn)), ...
                'FiniteDifferenceType','forward','ScaleProblem','jacobian');
            [yp,~,rp,efp,op]=lsqnonlin(@(v)fun(v),y,lower,upper,polish_opts);
            fn=fn+op.funcCount;
            % The exit flag describes optimizer termination, not the
            % constitutive contract.  Preserve any finite strict decrease;
            % the unchanged final residual assertion remains authoritative.
            if all(isfinite(rp)) && norm(rp,inf)<norm(residual,inf)
                y=yp; residual=rp; exitflag=efp;
            end
        catch
            % Keep the strict Newton result and let the final assertion
            % report the actual constitutive residual if polishing fails.
        end
    end
else
    options=optimoptions('lsqnonlin','Display','off','FunctionTolerance',1e-14, ...
        'StepTolerance',1e-14,'OptimalityTolerance',1e-13, ...
        'MaxIterations',get_option(closure,'phase_max_iterations',300), ...
        'MaxFunctionEvaluations',get_option(closure,'phase_max_function_evaluations',12000), ...
        'FiniteDifferenceType',get_option(closure,'finite_difference_type','central'), ...
        'ScaleProblem','jacobian');
    objective_scale=1e-6;
    [y,~,residual,exitflag,output]=lsqnonlin(@(v)fun(v)/objective_scale, ...
        y,lower,upper,options);
end
if transformed
    if strcmp(transform_kind,'logq')
        q_trial_final=[y(1);exp(y(2:end))];
        y=min(max(map_logq_to_y(y,sgn,scale_dt,mcoord), ...
            lower_phys),upper_phys);
    else
        q_trial_final=y;
        y=min(max(map_q_to_y(y,sgn,scale_dt,mcoord), ...
            lower_phys),upper_phys);
    end
    lower=lower_phys; upper=upper_phys;
    fun=fun_phys;
end
if ~get_option(closure,'strict_source_kinetics',false), q_trial_final=[]; end
[R,phase,terms]=evaluate(y,p,old,macro_strain,T,dt,closure,Dxx,xi0, ...
    xi_res,xi_recent,Btr,deta_p,Ybase,Ymult,resistance,capacity,direction_R,direction_M,w, ...
    population_mode,np,nx,pop_weights,kinetic_exponent_R,kinetic_exponent_M,kinetic_reference_rate, ...
    source_force_representation,Omega_state,q_trial_final);
rnorm=norm(R,inf);
% The natural-q coordinate is usually fastest, but near a projected
% activation/deactivation kink its semismooth Jacobian can lose the
% branch needed by the thermal trial.  This is a solver-coordinate issue,
% not a constitutive alternative: retry the identical equations in the
% logarithmic-q coordinate before rejecting the trial.  The fallback is
% fail-closed (it is only accepted after the original residual is rebuilt
% and passes the same strict tolerance).
if rnorm>get_option(closure,'residual_tolerance',2e-8) && ...
        strcmpi(solver_mode,'newton_q') && ...
        ~get_option(closure,'disable_logq_fallback',false)
    closure_fallback=closure;
    closure_fallback.phase_solver_mode='newton_logq';
    closure_fallback.disable_logq_fallback=true;
    % The logarithmic coordinate has a strictly positive realmin floor.  If
    % the source KKT branch is inactive, carrying the primary q trial into
    % that fallback can leave a finite pseudo-increment at the physical
    % lower bound (the m=200 onset case is especially sensitive).  Rebuild
    % the fallback from the exact admissible bound for inactive channels;
    % this is an active-set transaction, not an endpoint assignment.
    fallback_guess=y;
    for ii=1:n
        rt_active_for_direction=terms.kinetic_audit{ii}.source_condition_active && ...
            terms.kinetic_audit{ii}.source_branch==direction_R(ii);
        if ~rt_active_for_direction
            if direction_R(ii)>0, fallback_guess(1+ii)=lower_phys(1+ii);
            else, fallback_guess(1+ii)=upper_phys(1+ii); end
        end
        mt_active_for_direction=terms.kinetic_audit{ii}.MT.source_condition_active && ...
            terms.kinetic_audit{ii}.MT.source_branch==direction_M(ii);
        if ~mt_active_for_direction
            if direction_M(ii)>0, fallback_guess(1+n+ii)=lower_phys(1+n+ii);
            else, fallback_guess(1+n+ii)=upper_phys(1+n+ii); end
        end
    end
    % A failed natural-q iterate is not necessarily a valid seed in log-q:
    % the q coordinate can stop on a projected kink while its physical
    % increment is still several orders away from the logarithmic branch.
    % First rebuild the logarithmic solve from its own equation-derived
    % initial seed.  If that seed fails, retain the previous active-set
    % continuation as a second, fail-closed numerical attempt.  Neither
    % path changes the constitutive equations or accepts a non-strict
    % residual; both recursive calls must pass the same final gate.
    try
        [phase,audit]=sma.candidate_gradient_monolithic_phase_step_1d( ...
            p,old,macro_strain,T,dt,closure_fallback,[]);
    catch fallback_seed_exception
        try
            [phase,audit]=sma.candidate_gradient_monolithic_phase_step_1d( ...
                p,old,macro_strain,T,dt,closure_fallback,fallback_guess);
        catch fallback_guess_exception
            % Preserve the more informative second failure while retaining
            % the first identifier in the diagnostic message.
            error(fallback_guess_exception.identifier, ...
                '%s (natural log-q seed first failed: %s)', ...
                fallback_guess_exception.message,fallback_seed_exception.message);
        end
    end
    audit.solver_fallback='newton_q_to_newton_logq';
    audit.solver_primary_mode='newton_q';
    return;
end
if rnorm>get_option(closure,'residual_tolerance',2e-8)
    error('sma:candidate_gradient_monolithic_phase_step_1d:Closure', ...
        ['Monolithic gradient phase solve failed; residual %.6g, ' ...
        'exitflag %d, components %s.'],rnorm,exitflag,mat2str(R(:).',6));
end
audit=struct('stress_MPa',y(1),'delta_lambda_R',y(2:1+n), ...
    'delta_lambda_M',y(2+n:end),'residual',R,'residual_norm',rnorm, ...
    'terms',terms,'iterations',output.iterations,'function_count',output.funcCount, ...
    'optimizer_exitflag',exitflag, ...
    'operator',operator,'observation_accessed',false,'production_connected',false);
audit.phase_solver_mode=solver_mode;
audit.residual_event_lower_enabled=has_event_lower;
audit.residual_event_lower_applicable=has_event_lower & direction_M<0;
audit.residual_event_xi_M_lower=event_xi_lower;
audit.residual_event_boundary_active=has_event_lower && ...
    any(direction_M<0 & abs(y(2+n:end)- ...
    (event_xi_lower-xi0))<=2e-8);
if isfield(closure,'phase_debug') && closure.phase_debug
    at=y; hact=1e-7*max(1,abs(y));
    for kk=2:numel(y)
        if y(kk)-lower(kk)<=hact(kk), at(kk)=lower(kk); end
        if upper(kk)-y(kk)<=hact(kk), at(kk)=upper(kk); end
    end
    audit.debug_active_trial=at;
    audit.debug_active_residual=fun(at);
end
end

function [R,phase,terms]=evaluate(y,p,old,target,T,dt,c,Dxx,xi0, ...
    xi_res,xi_recent,Btr,deta_p,Ybase,Ymult,resistance,capacity,sR,sM,w, ...
    population_mode,np,nx,pop_weights,kinetic_exponent_R,kinetic_exponent_M,kinetic_reference_rate, ...
    source_force_representation,Omega_state,q_trial_override)
if nargin<30 || isempty(q_trial_override)
    q_trial_override=[];
end
n=numel(old); sigma=y(1); dR=y(2:1+n); dM=y(2+n:end);
qR_override=[]; qM_override=[];
if ~isempty(q_trial_override)
    assert(numel(q_trial_override)==1+2*n && all(isfinite(q_trial_override)), ...
        'The transformed kinetic override has an invalid width.');
    qR_override=max(q_trial_override(2:1+n),0);
    qM_override=max(q_trial_override(2+n:end),0);
end
phase=cell(n,1); inc=cell(n,1); strain=zeros(n,1);
piR=zeros(n,1); piM=zeros(n,1); rateR=zeros(n,1); rateM=zeros(n,1);
deta_mix=zeros(n,1);
path_strain=zeros(n,1); path_increment=zeros(n,1);
path_dissipation=zeros(n,1); path_energy_defect=zeros(n,1);
path_audit=cell(n,1);
thermal_strain=zeros(n,1); thermal_unlock=zeros(n,1);
thermal_dissipation=zeros(n,1); thermal_energy_defect=zeros(n,1);
thermal_audit=cell(n,1);
force_ledger=cell(n,1);
kinetic_audit=cell(n,1);
for i=1:n
    if strcmp(c.phase_composition_mode,'source_current_mixture')
        [phase{i},inc{i}]= ...
            sma.candidate_zhang_current_mixture_phase_increment_1d( ...
            p,old{i},dR(i),dM(i));
    else
        [phase{i},inc{i}]= ...
            sma.candidate_coupled_phase_increment(p,old{i},dR(i),dM(i));
    end
    mixed=max(phase{i}.xi_A+phase{i}.xi_R,eps);
    relR=phase{i}.xi_R/mixed;
    tauR=p.R_tau0+p.R_hardening*relR;
    piR(i)=p.R_H*sigma+p.R_entropy*(T(i)-c.Tref_K)-tauR;
    path_increment(i)=abs(dR(i))+abs(dM(i));
    if c.enable_path_accommodation
        [path_state,path_audit{i}]= ...
            sma.candidate_transformation_path_accommodation_1d( ...
            old{i}.path_accommodation_strain,sigma,path_increment(i), ...
            c.path_accommodation_moduli_MPa, ...
            c.path_accommodation_resistances_MPa);
        phase{i}.path_accommodation_strain=path_state;
        path_strain(i)=sum(path_state);
        path_dissipation(i)=path_audit{i}.dissipation_MPa;
        path_energy_defect(i)=path_audit{i}.energy_identity_defect_MPa;
    else
        path_audit{i}=struct('effective_strain',0, ...
            'dissipation_MPa',0,'energy_identity_defect_MPa',0, ...
            'path_increment',0,'observation_accessed',false, ...
            'production_connected',false);
    end
    if c.enable_thermal_accommodation
        % Use the RT KKT driving force, not dR>0.  The log-q solver
        % represents an inactive physical increment by a positive realmin
        % floor, so a raw increment test would create a false thermal path.
        % The default continuous weight also removes the finite jump of a
        % hard activity indicator at pi_R=0.
        gate_mode=get_option(c,'thermal_accommodation_gate','continuous_kkt');
        available=(sR(i)>0) && (old{i}.xi_A>100*eps(1));
        if strcmpi(gate_mode,'hard_kkt')
            rt_activity_weight=double(available && piR(i)>0);
        elseif strcmpi(gate_mode,'continuous_kkt')
            rt_activity_weight=double(available)*min(1,max(0,piR(i)/c.Y_R_MPa));
        else
            error('Unknown thermal-accommodation gate mode %s.',gate_mode);
        end
        [thermal_state,thermal_audit{i}]= ...
            sma.candidate_thermal_accommodation_1d( ...
            old{i}.thermal_accommodation_strain,sigma,rt_activity_weight, ...
            c.committed_temperature_K(i),T(i),p.R_entropy,c.Y_R_MPa, ...
            c.thermal_accommodation_moduli_MPa, ...
            c.thermal_accommodation_resistances_MPa);
        phase{i}.thermal_accommodation_strain=thermal_state;
        thermal_strain(i)=sum(thermal_state);
        thermal_unlock(i)=thermal_audit{i}.unlock_path;
        thermal_dissipation(i)=thermal_audit{i}.dissipation_MPa;
        thermal_energy_defect(i)=thermal_audit{i}.energy_identity_defect_MPa;
    else
        thermal_audit{i}=struct('effective_strain',0,'unlock_path',0, ...
            'dissipation_MPa',0,'energy_identity_defect_MPa',0, ...
            'forward_rt_active',false,'observation_accessed',false, ...
            'production_connected',false);
    end
    E=sma.effective_modulus(p,phase{i}.xi_M);
    strain(i)=sigma/E+p.alphaA*(T(i)-p.T0)+phase{i}.eps_RT+ ...
        phase{i}.eps_MT+path_strain(i)+thermal_strain(i);
end
xi=cellfun(@(s)s.xi_M,phase);
if population_mode && c.kappa_MPa_length2~=0
    [gradient,population_gradient_audit]= ...
        sma.candidate_population_gradient_map_1d( ...
        xi0,xi,np,nx,pop_weights,c.kappa_MPa_length2,Dxx);
elseif c.kappa_MPa_length2==0
    gradient=zeros(n,1);
    population_gradient_audit=struct('layout','zero-gradient', ...
        'observation_accessed',false,'production_connected',false);
else
    gradient=c.kappa_MPa_length2*Dxx*(0.5*(xi0+xi));
    population_gradient_audit=struct('layout','single-population', ...
        'observation_accessed',false,'production_connected',false);
end
for i=1:n
    % Keep the source-resistance implementation in one place.  The
    % gradient candidate is a Zhang-source shadow path, so its intrinsic
    % MT resistance must be identical to the audited Eqs. (61)--(63)
    % adapter rather than a second, potentially drifting formula.
    tauM=sma.candidate_zhang_source_mt_resistance( ...
        p,sM(i),phase{i}.xi_M,xi_recent(i),xi_res(i));
    % Zhang Eq.(37-b): eta_M minus the A/R parent-mixture entropy.
    deta=sma.mixed_mt_entropy(p,inc{i}.weight_R,deta_p(i));
    deta_mix(i)=deta;
    compliance=1/p.EM-1/p.EA;
    source_force=sma.candidate_zhang_source_force_representation_1d( ...
        inc{i}.g_mixed,sigma,Btr(i),Omega_state(i),source_force_representation);
    piM(i)=source_force.mechanical_force_MPa+0.5*sigma^2*compliance+ ...
        deta*(T(i)-c.Tref_K)-tauM+gradient(i);
    if isempty(q_trial_override)
        rateR(i)=sR(i)*c.Y_R_MPa*(abs(dR(i))/(dt*kinetic_reference_rate))^(1/kinetic_exponent_R);
        rateM(i)=sM(i)*(Ybase(i)*Ymult(i)+resistance(i))*(abs(dM(i))/(dt*kinetic_reference_rate))^(1/kinetic_exponent_M);
    else
        % In transformed coordinates the physical increment may underflow
        % even though q is finite (the first MT activation at m=200 is the
        % concrete example).  Preserve q for the source-rate residual while
        % retaining the physical dLambda for state/strain updates.
        rateR(i)=sR(i)*c.Y_R_MPa*qR_override(i);
        rateM(i)=sM(i)*(Ybase(i)*Ymult(i)+resistance(i))*qM_override(i);
    end
    elastic_force=0.5*sigma^2*compliance;
    thermal_force=deta*(T(i)-c.Tref_K);
    transformation_force=source_force.stress_transformation_MPa;
    reaction_force=source_force.main_Btr_contraction_MPa;
    selected_history_force=source_force.selected_history_force_MPa;
    known_force=transformation_force+selected_history_force+elastic_force+ ...
        thermal_force-tauM+gradient(i);
    force_ledger{i}=struct('elastic_MPa',elastic_force, ...
        'thermal_MPa',thermal_force, ...
        'transformation_MPa',transformation_force, ...
        'reaction_backstress_MPa',reaction_force, ...
        'g_mixed',inc{i}.g_mixed,'Btr_MPa',Btr(i), ...
        'reaction_backstress_reconstruction_MPa',inc{i}.g_mixed*Btr(i), ...
        'reaction_backstress_reconstruction_residual_MPa',reaction_force- ...
            inc{i}.g_mixed*Btr(i), ...
        'Btr_counting_policy','Eq36b_first_contraction_once', ...
        'reaction_backstress_selected',source_force.btr_selected, ...
        'source_force_representation',source_force.representation, ...
        'source_force_formula_id',source_force.formula_id, ...
        'selected_history_force_MPa',selected_history_force, ...
        'supplement_Omega_MPa',source_force.supplement_Omega_MPa, ...
        'supplement_minus_main_MPa',source_force.supplement_minus_main_MPa, ...
        'double_count_forbidden',source_force.double_count_forbidden, ...
        'intrinsic_hardening_MPa',-tauM, ...
        'intrinsic_hardening_source','candidate_zhang_source_mt_resistance', ...
        'gradient_MPa',gradient(i), ...
        'rate_MPa',NaN,'configurational_MPa',NaN, ...
        'common_constraint_MPa',NaN,'omega_MPa',NaN, ...
        'beta_companion_MPa',NaN, ...
        'source_literal_scope',source_literal_scope(source_force.representation), ...
        'source_literal_sum_MPa',known_force-gradient(i), ...
        'source_literal_residual_to_kernel_pi_MPa',-gradient(i), ...
        'source_literal_complete',gradient(i)==0, ...
        'source_literal_extension_terms_required',false, ...
        'extension_complete',false,'complete',false, ...
        'unresolved_scope','optional_project_extension', ...
        'known_sum_MPa',known_force,'known_sum_residual_MPa',known_force-piM(i), ...
        'unresolved',{{'rate_MPa','configurational_MPa', ...
        'common_constraint_MPa','omega_MPa','beta_companion_MPa'}});
end
scale=max(p.reference_stress,1);
strict_source_kinetics=logical(get_option(c,'strict_source_kinetics',false));
if strict_source_kinetics
    % q/KKT is the acceptance block in this opt-in path; do not add a
    % force-normalized quantity to a phase-fraction increment.
    rawR=zeros(n,1); rawM=zeros(n,1);
else
    rawR=(piR-rateR)/scale; rawM=(piM-rateM)/scale;
end
projectedR=zeros(n,1); projectedM=zeros(n,1);
loR_bound=zeros(n,1); hiR_bound=zeros(n,1);
loM_bound=zeros(n,1); hiM_bound=zeros(n,1);
projected_increment_R=zeros(n,1); projected_increment_M=zeros(n,1);
projection_receipt=cell(n,1);
for i=1:n
    loR=-old{i}.xi_R; hiR=old{i}.xi_A;
    loM=-old{i}.xi_M+xi_res(i); hiM=min(old{i}.xi_A+old{i}.xi_R, ...
        capacity(i)-old{i}.xi_M);
    if isfield(c,'residual_event_xi_M_lower') && sM(i)<0
        loM=max(loM,c.residual_event_xi_M_lower(i)-old{i}.xi_M);
    end
    if sR(i)>0, loR=0; else, hiR=0; end
    if sM(i)>0, loM=0; else, hiM=0; end
    loR_bound(i)=loR; hiR_bound(i)=hiR;
    loM_bound(i)=loM; hiM_bound(i)=hiM;
    projected_increment_R(i)=min(max(dR(i)+rawR(i),loR),hiR);
    projected_increment_M(i)=min(max(dM(i)+rawM(i),loM),hiM);
    projectedR(i)=dR(i)-projected_increment_R(i);
    projectedM(i)=dM(i)-projected_increment_M(i);
    if strict_source_kinetics
        % The historical projected map can call a numerically tiny
        % positive increment "lower-active" even when the source q-demand
        % is finite.  In the opt-in gate that is not an admissible kinetic
        % root: retain the force/q residual unless the *upper* capacity is
        % genuinely saturated.  This is a KKT diagnostic only; the default
        % production map above is unchanged.
        projectedR(i)=strict_source_projected_residual( ...
            piR(i),c.Y_R_MPa,sR(i),dR(i),dt,kinetic_reference_rate, ...
            kinetic_exponent_R,loR,hiR,old{i}.xi_A>0,old{i}.xi_R>0,scale, ...
            qR_override_value(qR_override,i));
        projectedM(i)=strict_source_projected_residual( ...
            piM(i),Ybase(i)*Ymult(i)+resistance(i),sM(i),dM(i),dt, ...
            kinetic_reference_rate,kinetic_exponent_M,loM,hiM, ...
            old{i}.xi_A+old{i}.xi_R>0,old{i}.xi_M>xi_res(i),scale, ...
            qR_override_value(qM_override,i));
    end
    projection_receipt{i}=struct( ...
        'contract','candidate_phase_projected_map_v1', ...
        'direction_R',sR(i),'direction_M',sM(i), ...
        'delta_R',dR(i),'delta_M',dM(i), ...
        'force_residual_R_MPa',piR(i)-rateR(i), ...
        'force_residual_M_MPa',piM(i)-rateM(i), ...
        'reference_stress_MPa',scale, ...
        'capacity_M',capacity(i),'xi_res',xi_res(i), ...
        'raw_increment_R',rawR(i),'raw_increment_M',rawM(i), ...
        'loR',loR,'hiR',hiR,'loM',loM,'hiM',hiM, ...
        'projected_increment_R',projected_increment_R(i), ...
        'projected_increment_M',projected_increment_M(i), ...
        'projected_residual_R',projectedR(i), ...
        'projected_residual_M',projectedM(i), ...
        'R_lower_active',abs(dR(i)-loR)<=100*eps(max(1,abs(loR))), ...
        'R_upper_active',abs(dR(i)-hiR)<=100*eps(max(1,abs(hiR))), ...
        'M_lower_active',abs(dM(i)-loM)<=100*eps(max(1,abs(loM))), ...
        'M_upper_active',abs(dM(i)-hiM)<=100*eps(max(1,abs(hiM))), ...
        'active_set_signature',sprintf('Rlo=%d|Rhi=%d|Mlo=%d|Mhi=%d', ...
            abs(dR(i)-loR)<=100*eps(max(1,abs(loR))), ...
            abs(dR(i)-hiR)<=100*eps(max(1,abs(hiR))), ...
            abs(dM(i)-loM)<=100*eps(max(1,abs(loM))), ...
            abs(dM(i)-hiM)<=100*eps(max(1,abs(hiM)))), ...
        'source_common_rt_complete',false, ...
        'observation_accessed',false,'production_connected',false);
    kinetic_audit{i}=kinetic_audit_receipt( ...
        piR(i),c.Y_R_MPa,sR(i),dR(i),dt,kinetic_reference_rate, ...
        kinetic_exponent_R,loR,hiR,rateR(i), ...
        projection_receipt{i}.R_lower_active,projection_receipt{i}.R_upper_active, ...
        old{i}.xi_A>0,old{i}.xi_R>0,'RT',c, ...
        qR_override_value(qR_override,i)); %#ok<AGROW>
    kinetic_audit{i}.MT=kinetic_audit_receipt( ...
        piM(i),Ybase(i)*Ymult(i)+resistance(i),sM(i),dM(i),dt, ...
        kinetic_reference_rate,kinetic_exponent_M,loM,hiM,rateM(i), ...
        projection_receipt{i}.M_lower_active,projection_receipt{i}.M_upper_active, ...
        old{i}.xi_A+old{i}.xi_R>0,old{i}.xi_M>xi_res(i),'MT',c, ...
        qR_override_value(qM_override,i)); %#ok<AGROW>
end
control_residual=(w.'*strain-target)/max(p.Hmax_i,1e-4);
if strcmp(get_option(c,'candidate_mechanical_control','prescribed_strain'), ...
        'validation_prescribed_stress')
    % One MPa is the explicit stress residual normalization. Final stress
    % control is independently checked against the unchanged 2e-7 MPa gate.
    control_residual=sigma-c.candidate_target_stress_MPa;
end
R=[control_residual;projectedR;projectedM];
terms=struct('cell_strain',strain,'macro_strain',w.'*strain, ...
    'gradient_driving_MPa',gradient, ...
    'population_gradient_audit',population_gradient_audit, ...
    'pi_R_MPa',piR,'pi_M_MPa',piM, ...
    'rate_R_MPa',rateR,'rate_M_MPa',rateM,'xi_M',xi);
terms.force_ledger=force_ledger;
terms.phase_projection_receipt=projection_receipt;
terms.raw_increment_R=rawR; terms.raw_increment_M=rawM;
terms.kinetic_audit=kinetic_audit;
terms.loR=loR_bound; terms.hiR=hiR_bound;
terms.loM=loM_bound; terms.hiM=hiM_bound;
terms.projected_increment_R=projected_increment_R;
terms.projected_increment_M=projected_increment_M;
terms.projected_residual_R=projectedR;
terms.projected_residual_M=projectedM;
receipt_required=get_option(c,'require_zhang_source_receipt',false);
terms.source_receipt_required=receipt_required;
terms.source_receipt=conditional_source_receipt(c,receipt_required);
history_receipt_required_eval=get_option(c,'require_history_state_receipt',false);
terms.history_state_receipt_required=history_receipt_required_eval;
terms.history_state_receipt=conditional_history_receipt(c,history_receipt_required_eval);
terms.path_accommodation_active=logical(c.enable_path_accommodation);
terms.path_accommodation_strain=path_strain;
terms.path_accommodation_increment=path_increment;
terms.path_accommodation_dissipation_MPa=path_dissipation;
terms.path_accommodation_energy_defect_MPa=path_energy_defect;
terms.path_accommodation_audit=path_audit;
terms.thermal_accommodation_active=logical(c.enable_thermal_accommodation);
terms.thermal_accommodation_strain=thermal_strain;
terms.thermal_accommodation_unlock=thermal_unlock;
terms.thermal_accommodation_dissipation_MPa=thermal_dissipation;
terms.thermal_accommodation_energy_defect_MPa=thermal_energy_defect;
terms.thermal_accommodation_audit=thermal_audit;
terms.delta_eta_mix_MPaK=deta_mix;
terms.resistance_multiplier=Ymult;
terms.mt_resistance_source='candidate_zhang_source_mt_resistance';
terms.mt_kinetic_threshold_source='Yf0_plus_d_xi_recent_power_minus_Yp';
terms.transformation_exponent_R=kinetic_exponent_R;
terms.transformation_exponent_M=kinetic_exponent_M;
terms.transformation_reference_rate_per_s=kinetic_reference_rate;
terms.strict_source_kinetics=strict_source_kinetics;
if strict_source_kinetics
    terms.acceptance_residual_kind='strain_over_H_and_dimensionless_projected_q';
else
    terms.acceptance_residual_kind='strain_over_H_and_historical_projected_increment';
end
terms.phase_composition_mode=c.phase_composition_mode;
end

function scope=source_literal_scope(representation)
switch representation
    case 'main_eq36b'
        scope='Zhang_Eq36b';
    case 'supplement_s4a'
        scope='Zhang_Supplement_S4a';
    otherwise
        error('Unknown Zhang source-force representation: %s.',representation);
end
end

function value=conditional_source_receipt(closure,required)
if required
    value=closure.Y_M_threshold_receipt;
else
    value={};
end
end

function a=canonical_eq39_ablation()
a=struct('enabled',false,'m_M',200,'m_R',200,'reference_rate_per_s',1, ...
    'ablation_id','canonical','parent_source_receipt_digest','canonical', ...
    'source_authoritative',true,'production_connected',false, ...
    'experimental_acceptance',false,'observation_accessed',false, ...
    'experimental_backfill',false,'tag_suffix','canonical');
end

function value=conditional_history_receipt(closure,required)
if required
    value=closure.history_state_receipt;
else
    value={};
end
end

function token=history_state_token(receipt)
% Deterministic identity for the exact history packet consumed by this
% kernel.  It is intentionally a readable token, not a cryptographic hash.
ablation_id='canonical'; exponent=200;
if isfield(receipt,'eq39_ablation_id'), ablation_id=receipt.eq39_ablation_id; end
if isfield(receipt,'eq39_exponent'), exponent=receipt.eq39_exponent; end
token=sprintf(['gradient_history_v2|case=%d|xi_Mr=%.17g|Omega=%.17g|Yp=%.17g|' ...
    'deta_p=%.17g|Btr=%.17g|Yf=%.17g|Ybase=%.17g|card=%s|formula=%s'], ...
    receipt.material_case_id,receipt.xi_Mr,receipt.Omega_MPa,receipt.Yp_MPa, ...
    receipt.DeltaEtaP_MPaK,receipt.Btr_MPa,receipt.Yf_MPa,receipt.Y_M_base_MPa, ...
    receipt.source_card_id,receipt.formula_id);
token=sprintf('%s|eq39_m=%.17g|ablation=%s',token,exponent,ablation_id);
if isfield(receipt,'source_semantic_card_id')
    token=sprintf('%s|branch=%s',token,receipt.source_semantic_card_id);
end
end

function value=expand(value,n)
if isscalar(value), value=repmat(value,n,1); else, value=value(:); end
assert(numel(value)==n && all(isfinite(value)),'One finite value per cell is required.');
end

function value=get_option(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end

function receipt=kinetic_audit_receipt(pi,y,dir,d,dt,r0,m,lo,hi,rate, ...
    lower_active,upper_active,available_forward,available_reverse,channel,closure, ...
    q_trial_override)
%KINETIC_AUDIT_RECEIPT Observation-free source-rate diagnostics.
% The fields below do not enter the residual or state update.  They expose
% the distinction between the printed rate law in q-space and the absolute
% phase-fraction projected residual, which is essential when m=200.
assert(y>0 && dt>0 && r0>0 && m>0,'Invalid kinetic audit scale.');
scale=dt*r0;
if nargin<17 || isempty(q_trial_override)
    q_trial=0;
    if abs(d)>0, q_trial=(abs(d)/scale)^(1/m); end
else
    q_trial=max(q_trial_override,0);
end
source_sign=0;
if available_forward && pi>0, source_sign=1;
elseif available_reverse && pi<0, source_sign=-1;
end
source_active=source_sign~=0;
q_demand=0;
if source_active, q_demand=abs(pi)/y; end
source_rate_force=source_sign*y*q_demand;
source_increment=source_sign*scale*q_demand^m;
% Use strict physical inequalities here.  The projected residual's
% 100*eps activity band is intentionally not reused as a material zero.
physical_interior=(d>lo && d<hi);
q_residual=NaN;
if source_active && (physical_interior || ~isempty(q_trial_override))
    q_residual=q_trial-q_demand;
end
q_upper=(max(abs([lo,hi]))/scale)^(1/m);
coordinate_kind=char(get_option(closure,'phase_solver_mode','direct'));
if strcmpi(coordinate_kind,'newton_logq')
    q_floor=exp(log(realmin('double'))/m);
    mapped_floor_increment=scale*realmin('double');
elseif strcmpi(coordinate_kind,'newton_q')
    q_floor=0;
    mapped_floor_increment=0;
else
    q_floor=0;
    mapped_floor_increment=0;
end
receipt=struct( ...
    'channel',channel,'coordinate_kind',coordinate_kind, ...
    'dt_s',dt,'reference_rate_per_s',r0,'exponent',m, ...
    'pi_MPa',pi,'threshold_MPa',y,'direction',dir, ...
    'source_phase_available_forward',logical(available_forward), ...
    'source_phase_available_reverse',logical(available_reverse), ...
    'source_branch',source_sign,'source_condition_active',logical(source_active), ...
    'q_demand',q_demand,'q_trial',q_trial, ...
    'q_bound_lower',0,'q_bound_upper',q_upper, ...
    'q_residual',q_residual, ...
    'source_rate_force_MPa',source_rate_force, ...
    'trial_rate_force_MPa',rate, ...
    'trial_rate_force_from_q_MPa',dir*y*q_trial, ...
    'source_increment',source_increment, ...
    'trial_increment',d,'increment_residual',d-source_increment, ...
    'physical_lower_bound',lo,'physical_upper_bound',hi, ...
    'physical_interior',logical(physical_interior), ...
    'physical_zero_tolerance',0,'underflow_to_zero',logical(source_increment==0 && source_active), ...
    'lower_bound_active_by_projection',logical(lower_active), ...
    'upper_bound_active_by_projection',logical(upper_active), ...
    'physical_active_set_signature',sprintf('%s|src=%d|int=%d|lo=%d|hi=%d', ...
        channel,source_sign,physical_interior,lower_active,upper_active), ...
    'q_floor',q_floor,'mapped_floor_increment',mapped_floor_increment, ...
    'rate_force_residual_MPa',pi-dir*y*q_trial, ...
    'source_rate_residual_MPa',pi-source_rate_force, ...
    'source_exact_interior_residual_MPa',q_residual, ...
    'observation_accessed',false,'production_connected',false);
end

function residual=strict_source_projected_residual(pi,y,dir,d,dt,r0,m, ...
    lo,hi,available_forward,available_reverse,scale,q_trial_override) %#ok<INUSD>
%STRICT_SOURCE_PROJECTED_RESIDUAL Separate dimensionless q/KKT block.
if nargin>=13 && ~isempty(q_trial_override)
    q_trial=max(q_trial_override,0);
else
    q_trial=0;
    if abs(d)>0
        q_trial=(abs(d)/(dt*r0))^(1/m);
    end
end
available=(dir>0 && available_forward) || (dir<0 && available_reverse);
residual=sma.candidate_source_q_kkt_residual_1d( ...
    pi,y,dir,q_trial,dt,r0,m,max(abs([lo,hi])),available);
end

function z=map_y_to_q(y,sgn,scale_dt,m)
z=y; z(1)=y(1);
for k=2:numel(y)
    z(k)=(abs(y(k))/scale_dt)^(1/m(k-1));
end
end

function R=evaluate_q_coordinate(v,kind,sgn,scale_dt,mcoord, ...
    lower_phys,upper_phys,eval_args)
%EVALUATE_Q_COORDINATE Preserve q through an underflowing dLambda map.
if strcmpi(kind,'logq')
    q_trial=[v(1);exp(v(2:end))];
    physical=map_logq_to_y(v,sgn,scale_dt,mcoord);
else
    q_trial=v;
    physical=map_q_to_y(v,sgn,scale_dt,mcoord);
end
physical=min(max(physical,lower_phys),upper_phys);
[R,~,~]=evaluate(physical,eval_args{:},q_trial);
end

function y=map_q_to_y(z,sgn,scale_dt,m)
y=z; y(1)=z(1);
for k=2:numel(z)
    y(k)=sgn(k-1)*scale_dt*max(z(k),0)^m(k-1);
end
end

function z=map_y_to_logq(y,sgn,scale_dt,m) %#ok<INUSD>
z=y; z(1)=y(1);
floor_value=realmin('double');
for k=2:numel(y)
    z(k)=log(max(abs(y(k))/scale_dt,floor_value))/m(k-1);
end
end

function y=map_logq_to_y(z,sgn,scale_dt,m)
y=z; y(1)=z(1);
for k=2:numel(z)
    y(k)=sgn(k-1)*scale_dt*exp(m(k-1)*z(k));
end
end

function value=qR_override_value(q,k)
if isempty(q)
    value=[];
else
    value=q(k);
end
end
