function [phase,plastic,T,audit]= ...
    candidate_gradient_monolithic_stress_thermomechanical_step_1d( ...
    p,old_phase,old_plastic,T_old,sigma_old,target_stress,dt,closure, ...
    plastic_closure,thermal_closure,old_macro_strain,guess)
%SMA.CANDIDATE_GRADIENT_MONOLITHIC_STRESS_THERMOMECHANICAL_STEP_1D
% Stress-control wrapper around the accepted gradient thermo-mechanical
% transaction. The outer scalar root adjusts macro strain; every probe
% solves temperature and phase/plastic state from the same committed state.

assert(isscalar(target_stress) && isfinite(target_stress) && ...
    isscalar(sigma_old) && isfinite(sigma_old) && isfinite(dt) && dt>0, ...
    'sma:candidate_gradient_monolithic_stress_thermomechanical_step_1d:Control', ...
    'Stress controls and time step must be finite scalars.');
assert(isscalar(old_macro_strain) && isfinite(old_macro_strain), ...
    'sma:candidate_gradient_monolithic_stress_thermomechanical_step_1d:OldStrain', ...
    'The previously accepted macro strain is required.');
if nargin<12 || isempty(guess), guess=struct(); end
if ~isstruct(guess), error('sma:candidate_gradient_monolithic_stress_thermomechanical_step_1d:Guess', ...
        'Guess must be a struct.'); end
tol=get_option(closure,'stress_control_tolerance_MPa',2e-7);
max_expand=get_option(closure,'stress_control_max_bracket_expansions',16);
max_iter=get_option(closure,'stress_control_max_iterations',60);
width_cap=get_option(closure,'stress_control_initial_width',0.005);
width_policy=char(get_option(closure,'stress_control_initial_width_policy','fixed'));
if strcmpi(width_policy,'elastic_increment')
    w=expand_local(closure.weights,numel(old_phase));
    Ecell=cellfun(@(state)sma.effective_modulus(p,state.xi_M),old_phase(:));
    Eeffective=1/sum(w./Ecell);
    width_floor=get_option(closure,'stress_control_initial_width_floor',1e-5);
    width_factor=get_option(closure,'stress_control_elastic_width_factor',3);
    width=min(width_cap,max(width_floor,width_factor* ...
        abs(target_stress-sigma_old)/Eeffective));
elseif strcmpi(width_policy,'fixed')
    width=width_cap;
else
    error('Unknown stress-control initial-width policy %s.',width_policy);
end
initial_width=width;
bracket_policy=char(get_option(closure,'stress_control_bracket_policy','symmetric'));
thermal_tol=get_option(thermal_closure,'joint_temperature_tolerance_K',1e-7);
require_signature=get_option(closure,'stress_control_require_single_signature',false);
probe_continuation=get_option(closure,'stress_control_probe_continuation',false);
residual_event_policy=char(get_option(plastic_closure,'residual_m_event_policy','clip_legacy'));
residual_event_rejection_mode=strcmpi(residual_event_policy,'reject_crossing');
residual_event_outer_policy=char(get_option(closure,'residual_event_outer_policy','throw'));
assert(any(strcmpi(residual_event_outer_policy,{'throw','return_boundary'})), ...
    'sma:candidate_gradient_monolithic_stress_thermomechanical_step_1d:ResidualEventOuterPolicy', ...
    'Unknown residual-event outer policy: %s.',residual_event_outer_policy);
assert(isfinite(tol) && tol>0 && max_expand>=1 && max_iter>=4 && ...
    isfinite(width) && width>0 && isfinite(thermal_tol) && thermal_tol>0, ...
    'sma:candidate_gradient_monolithic_stress_thermomechanical_step_1d:Closure', ...
    'Stress/thermal solver options are inadmissible.');
baseline_failure_identifier=''; baseline_failure_message='';
continuation_used=false;
if probe_continuation
    % Preserve the historical root path whenever it already satisfies all
    % gates.  Branch continuation is a fail-closed numerical fallback only
    % for outer stress-root failures; it must not hide constitutive/source
    % receipt errors raised below that layer.
    baseline_closure=closure;
    baseline_closure.stress_control_probe_continuation=false;
    try
        [phase,plastic,T,audit]= ...
            sma.candidate_gradient_monolithic_stress_thermomechanical_step_1d( ...
            p,old_phase,old_plastic,T_old,sigma_old,target_stress,dt, ...
            baseline_closure,plastic_closure,thermal_closure, ...
            old_macro_strain,guess);
        audit.probe_continuation_enabled=true;
        audit.probe_continuation_requested=true;
        audit.probe_continuation_used=false;
        audit.baseline_failure_identifier='';
        audit.baseline_failure_message='';
        return;
    catch baseline_exception
        recoverable={ ...
            'sma:candidate_gradient_monolithic_stress_thermomechanical_step_1d:Bracket', ...
            'sma:candidate_gradient_monolithic_stress_thermomechanical_step_1d:InvalidProbe', ...
            'sma:candidate_gradient_monolithic_stress_thermomechanical_step_1d:Convergence', ...
            'sma:candidate_gradient_monolithic_stress_thermomechanical_step_1d:Signature'};
        if ~any(strcmp(baseline_exception.identifier,recoverable))
            rethrow(baseline_exception);
        end
        baseline_failure_identifier=baseline_exception.identifier;
        baseline_failure_message=baseline_exception.message;
        continuation_used=true;
    end
end
centre=old_macro_strain;
if isfield(guess,'macro_strain') && isscalar(guess.macro_strain) && isfinite(guess.macro_strain)
    centre=guess.macro_strain;
end
probe_count=0; invalid_probe_count=0; continuation_seed_count=0;
last_invalid_identifier=''; last_invalid_message='';
last_residual_event_rejection=false; last_residual_event_message='';
[f0,~,~,~,a0]=probe(centre,guess);
if a0.valid_probe && abs(f0)<=tol
    phase=a0.phase_state; plastic=a0.plastic_state; T=a0.temperature_K;
    audit=assemble(a0,centre,f0,centre,centre,0,true,0);
    audit.probe_count=probe_count; audit.invalid_probe_count=invalid_probe_count;
    audit.initial_width=initial_width; audit.initial_width_policy=width_policy;
    audit.bracket_policy=bracket_policy;
    audit.probe_continuation_enabled=probe_continuation;
    audit.probe_continuation_requested=probe_continuation;
    audit.probe_continuation_used=continuation_used;
    audit.baseline_failure_identifier=baseline_failure_identifier;
    audit.baseline_failure_message=baseline_failure_message;
    audit.continuation_seed_count=continuation_seed_count;
    audit.last_invalid_identifier=last_invalid_identifier;
    audit.last_invalid_message=last_invalid_message;
    return;
end
expansions=0;
if strcmpi(bracket_policy,'directional_from_committed') && valid_probe(a0)
    if f0>0
        hi=centre; fhi=f0; ahi=a0;
        lo=centre-width; [flo,~,~,~,alo]=probe(lo,select_seed(lo,centre,centre,a0,a0));
        while ~(valid_probe(alo) && signed(flo)*signed(fhi)<=0) && ...
                expansions<max_expand
            previous_lo=lo; previous_alo=alo;
            width=2*width; lo=centre-width;
            [flo,~,~,~,alo]=probe(lo, ...
                select_seed(lo,previous_lo,centre,previous_alo,a0));
            expansions=expansions+1;
        end
    else
        lo=centre; flo=f0; alo=a0;
        hi=centre+width; [fhi,~,~,~,ahi]=probe(hi,select_seed(hi,centre,centre,a0,a0));
        while ~(valid_probe(ahi) && signed(flo)*signed(fhi)<=0) && ...
                expansions<max_expand
            previous_hi=hi; previous_ahi=ahi;
            width=2*width; hi=centre+width;
            [fhi,~,~,~,ahi]=probe(hi, ...
                select_seed(hi,centre,previous_hi,a0,previous_ahi));
            expansions=expansions+1;
        end
    end
elseif strcmpi(bracket_policy,'symmetric')
    lo=centre-width; hi=centre+width;
    [flo,~,~,~,alo]=probe(lo,select_seed(lo,centre,centre,a0,a0));
    [fhi,~,~,~,ahi]=probe(hi,select_seed(hi,centre,centre,a0,a0));
    while ~(valid_probe(alo) && valid_probe(ahi) && ...
            signed(flo)*signed(fhi)<=0) && expansions<max_expand
        previous_lo=lo; previous_hi=hi;
        previous_alo=alo; previous_ahi=ahi;
        width=2*width; lo=centre-width; hi=centre+width;
        [flo,~,~,~,alo]=probe(lo, ...
            select_seed(lo,previous_lo,previous_hi,previous_alo,previous_ahi));
        [fhi,~,~,~,ahi]=probe(hi, ...
            select_seed(hi,previous_lo,previous_hi,previous_alo,previous_ahi));
        expansions=expansions+1;
    end
else
    error('Unknown stress-control bracket policy %s.',bracket_policy);
end
if ~(valid_probe(alo) && valid_probe(ahi) && signed(flo)*signed(fhi)<=0)
    if residual_event_rejection_mode && last_residual_event_rejection
        if strcmpi(residual_event_outer_policy,'return_boundary')
            phase=old_phase; plastic=old_plastic; T=T_old;
            audit=struct('schema_id','NITI_S4_RESIDUAL_EVENT_BOUNDARY_AUDIT_V1', ...
                'valid_probe',false,'event_boundary',true, ...
                'residual_event_rejection_mode',true, ...
                'residual_event_rejection_seen',true, ...
                'event_boundary_identifier','sma:candidate_gradient_monolithic_stress_thermomechanical_step_1d:ResidualEventBoundary', ...
                'event_boundary_message',last_residual_event_message, ...
                'target_stress_MPa',target_stress,'stress_MPa',NaN, ...
                'parent_preserved',true,'state_mutated',false, ...
                'observation_accessed',false,'production_connected',false);
            return;
        end
        error('sma:candidate_gradient_monolithic_stress_thermomechanical_step_1d:ResidualEventBoundary', ...
            ['Strict residual-M trial rejection reached the outer stress root ', ...
            'without an event-aware bracket: %s'],last_residual_event_message);
    end
    error('sma:candidate_gradient_monolithic_stress_thermomechanical_step_1d:Bracket', ...
        ['Could not bracket target stress with two valid probes; last invalid ' ...
        'probe %s: %s'],last_invalid_identifier,last_invalid_message);
end
best=alo; best_x=lo; best_f=flo; iter=0;
for iter=1:max_iter
    denom=fhi-flo;
    if abs(denom)<=eps(max(abs([fhi,flo])))
        x=0.5*(lo+hi);
    else
        x=min(max((lo*fhi-hi*flo)/denom,lo),hi);
    end
    [fx,~,~,~,ax]=probe(x,select_seed(x,lo,hi,alo,ahi));
    if ~valid_probe(ax)
        % A rejected thermal/mechanical probe is never a bracket endpoint.
        % Try interior cutback points from the same committed bracket before
        % failing closed; no state or heat ledger is written by these probes.
        [fx,x,ax]=recover_valid_probe(lo,hi,alo,ahi);
        if ~valid_probe(ax)
            error('sma:candidate_gradient_monolithic_stress_thermomechanical_step_1d:InvalidProbe', ...
                'No valid stress probe remained after local bracket cutback.');
        end
    end
    if abs(fx)<abs(best_f), best=ax; best_x=x; best_f=fx; end
    if abs(fx)<=tol
        best=ax; best_x=x; best_f=fx; break;
    end
    if signed(flo)*signed(fx)<=0
        hi=x; fhi=fx; ahi=ax;
    else
        lo=x; flo=fx; alo=ax;
    end
end
% Regula-falsi can stagnate when the coupled thermal solve contributes a
% small, non-smooth numerical perturbation to the scalar stress response.
% Preserve the strict stress tolerance and the same committed-state probe
% transaction, but finish the already valid sign bracket with safeguarded
% bisection rather than accepting a relaxed residual.
if abs(best_f)>tol && valid_probe(alo) && valid_probe(ahi) && ...
        signed(flo)*signed(fhi)<=0
    for jj=1:max(20,ceil(log2(max(1,(hi-lo)/eps(max(1,abs(best_x)))))))
        if jj>64, break; end
        x=0.5*(lo+hi);
        [fx,~,~,~,ax]=probe(x,select_seed(x,lo,hi,alo,ahi));
        if ~valid_probe(ax), continue; end
        if abs(fx)<abs(best_f), best=ax; best_x=x; best_f=fx; end
        if abs(fx)<=tol, break; end
        if signed(flo)*signed(fx)<=0
            hi=x; fhi=fx; ahi=ax;
        else
            lo=x; flo=fx; alo=ax;
        end
    end
end
% At a projected phase boundary, repeated bracketing can reach machine-scale
% strain width while the nested thermal/phase solve still varies by a few
% stress ulps according to its initial branch.  When explicit probe
% continuation is enabled, finish with a local same-branch Newton polish.
% This changes only the trial initial guess; the original constitutive and
% stress residuals remain the acceptance gate.
if abs(best_f)>tol && probe_continuation && valid_probe(best)
    for jj=1:12
        h=max(1e-9,min(1e-7,0.05*max(initial_width,hi-lo)));
        [fL,~,~,~,aL]=probe(best_x-h,best.continuation_guess);
        [fU,~,~,~,aU]=probe(best_x+h,best.continuation_guess);
        if ~(valid_probe(aL) && valid_probe(aU)), break; end
        slope_local=(fU-fL)/(2*h);
        if ~isfinite(slope_local) || abs(slope_local)<=eps, break; end
        xtrial=best_x-best_f/slope_local;
        [ftrial,~,~,~,atrial]=probe(xtrial,best.continuation_guess);
        if ~valid_probe(atrial) || abs(ftrial)>=abs(best_f), break; end
        best=atrial; best_x=xtrial; best_f=ftrial;
        if abs(best_f)<=tol, break; end
    end
end
assert(abs(best_f)<=tol, ...
    'sma:candidate_gradient_monolithic_stress_thermomechanical_step_1d:Convergence', ...
    ['Stress residual %.9g MPa exceeds tolerance %.9g MPa at strain %.17g; ' ...
    'bracket [%.17g, %.17g], endpoint residuals [%.9g, %.9g], probes %d.'], ...
    best_f,tol,best_x,lo,hi,flo,fhi,probe_count);
if require_signature
    h=max(1e-8,1e-4*max(1,abs(best_x)));
    [fL,~,~,~,aL]=probe(best_x-h,best.continuation_guess);
    [fU,~,~,~,aU]=probe(best_x+h,best.continuation_guess);
    assert(valid_probe(aL) && valid_probe(aU) && ...
        same_signature(aL.signature,best.signature) && ...
        same_signature(aU.signature,best.signature), ...
        'sma:candidate_gradient_monolithic_stress_thermomechanical_step_1d:Signature', ...
        'Stress root crosses an invalid or changing phase/capacity signature.');
    slope=(fU-fL)/(2*h);
    assert(isfinite(slope) && abs(slope)>0, ...
        'sma:candidate_gradient_monolithic_stress_thermomechanical_step_1d:Slope', ...
        'Stress root has no finite local monotonicity slope.');
else
    slope=NaN; aL=struct(); aU=struct();
end
phase=best.phase_state; plastic=best.plastic_state; T=best.temperature_K;
audit=assemble(best,best_x,best_f,lo,hi,iter,false,expansions);
audit.probe_count=probe_count; audit.invalid_probe_count=invalid_probe_count;
audit.initial_width=initial_width; audit.initial_width_policy=width_policy;
audit.bracket_policy=bracket_policy;
audit.probe_continuation_enabled=probe_continuation;
audit.probe_continuation_requested=probe_continuation;
audit.probe_continuation_used=continuation_used;
audit.baseline_failure_identifier=baseline_failure_identifier;
audit.baseline_failure_message=baseline_failure_message;
audit.continuation_seed_count=continuation_seed_count;
audit.last_invalid_identifier=last_invalid_identifier;
audit.last_invalid_message=last_invalid_message;
audit.residual_event_rejection_mode=residual_event_rejection_mode;
audit.residual_event_rejection_seen=last_residual_event_rejection;
audit.signature_gate_required=require_signature;
audit.local_stress_slope_MPa_per_strain=slope;
if require_signature
    audit.signature_left=aL.signature; audit.signature_right=aU.signature;
end

    function [residual,phase_state,plastic_state,temperature_K,probe_audit]=probe(macro_strain,local_guess)
        probe_count=probe_count+1;
        if nargin<2 || isempty(local_guess), local_guess=guess; end
        try
            [phase_state,plastic_state,temperature_K,full_audit]= ...
                sma.candidate_gradient_monolithic_joint_thermomechanical_step_1d( ...
                p,old_phase,old_plastic,T_old,sigma_old,macro_strain,dt,closure, ...
                plastic_closure,thermal_closure,local_guess);
            assert(isfinite(full_audit.stress_MPa) && ...
                full_audit.max_temperature_residual_K<=thermal_tol, ...
                'sma:candidate_gradient_monolithic_stress_thermomechanical_step_1d:Thermal', ...
                'Thermo-mechanical probe did not satisfy the thermal residual.');
            residual=full_audit.stress_MPa-target_stress;
            probe_audit=struct();
            probe_audit.phase_state=phase_state;
            probe_audit.plastic_state=plastic_state;
            probe_audit.temperature_K=temperature_K;
            probe_audit.full_audit=full_audit;
            probe_audit.stress_MPa=full_audit.stress_MPa;
            probe_audit.valid_probe=true;
            probe_audit.failure_identifier='';
            probe_audit.failure_message='';
            probe_audit.signature=make_signature(phase_state,full_audit,closure);
            probe_audit.continuation_guess=make_continuation_guess( ...
                temperature_K,full_audit);
        catch exception
            invalid_probe_count=invalid_probe_count+1;
            last_invalid_identifier=exception.identifier;
            last_invalid_message=exception.message;
            residual=NaN; phase_state=[]; plastic_state=[]; temperature_K=[];
            probe_audit=struct();
            probe_audit.phase_state=[]; probe_audit.plastic_state=[];
            probe_audit.temperature_K=[]; probe_audit.full_audit=struct( ...
                'observation_accessed',false,'production_connected',false);
            probe_audit.stress_MPa=NaN; probe_audit.valid_probe=false;
            probe_audit.failure_identifier=exception.identifier;
            probe_audit.failure_message=exception.message;
            probe_audit.residual_event_rejected=residual_event_rejection_mode && ...
                is_residual_event_exception(exception);
            if probe_audit.residual_event_rejected
                last_residual_event_rejection=true;
                last_residual_event_message=exception.message;
            end
            probe_audit.signature=struct('key','INVALID');
            probe_audit.continuation_guess=struct();
        end
    end

    function [residual,xout,aout]=recover_valid_probe(a,b,aa,ab)
        fractions=[0.5,0.25,0.75,0.125,0.875];
        residual=NaN; xout=NaN; aout=struct('valid_probe',false, ...
            'signature',struct('key','INVALID'));
        for jj=1:numel(fractions)
            xt=a+(b-a)*fractions(jj);
            [ft,~,~,~,at]=probe(xt,select_seed(xt,a,b,aa,ab));
            if valid_probe(at)
                residual=ft; xout=xt; aout=at; return;
            end
        end
    end

    function seed=select_seed(x,a,b,aa,ab)
        seed=guess;
        if ~probe_continuation, return; end
        if valid_probe(aa) && valid_probe(ab)
            if strcmpi(bracket_policy,'directional_from_committed')
                if target_stress<sigma_old, selected=ab; else, selected=aa; end
            elseif abs(x-a)<=abs(x-b)
                selected=aa;
            else
                selected=ab;
            end
        elseif valid_probe(aa)
            selected=aa;
        elseif valid_probe(ab)
            selected=ab;
        else
            return;
        end
        if isfield(selected,'continuation_guess') && ...
                ~isempty(fieldnames(selected.continuation_guess))
            seed=selected.continuation_guess;
            continuation_seed_count=continuation_seed_count+1;
        end
    end
end

function audit=assemble(probe_audit,macro_strain,residual,lo,hi,iterations,initial,expansions)
audit=struct();
audit.target_stress_MPa=probe_audit.stress_MPa-residual;
audit.stress_MPa=probe_audit.stress_MPa;
audit.stress_control_residual_MPa=residual;
audit.macro_strain=macro_strain;
audit.bracket_low=lo; audit.bracket_high=hi;
audit.iterations=iterations; audit.bracket_expansions=expansions;
audit.initial_probe_accepted=initial;
audit.full=probe_audit.full_audit;
audit.temperature_K=probe_audit.temperature_K;
audit.valid_probe=probe_audit.valid_probe;
audit.signature=probe_audit.signature;
audit.observation_accessed=false; audit.production_connected=false;
end

function value=get_option(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end

function value=signed(value)
if value>0, value=1; elseif value<0, value=-1; else, value=0; end
end

function tf=valid_probe(a)
tf=isstruct(a) && isfield(a,'valid_probe') && a.valid_probe && ...
    isfield(a,'stress_MPa') && isfinite(a.stress_MPa) && ...
    isfield(a,'temperature_K') && ~isempty(a.temperature_K) && ...
    isfield(a,'full_audit');
end

function tf=is_residual_event_exception(exception)
% Preserve the nested MException cause chain so the outer stress root can
% distinguish a strict residual-M moving-boundary rejection from an ordinary
% thermal/Newton invalid probe.
tf=false;
if ~isa(exception,'MException'), return; end
tf=contains(char(exception.identifier),'ResidualEventCrossing') || ...
    contains(char(exception.message),'ResidualEventCrossing');
if tf, return; end
causes=exception.cause;
for k=1:numel(causes)
    if is_residual_event_exception(causes{k}), tf=true; return; end
end
end

function sig=make_signature(phase_state,full_audit,closure)
n=numel(phase_state); xi=cellfun(@(s)s.xi_M,phase_state(:));
cap=expand_local(get_option(closure,'f_transformable',1),n);
tol=1e-12;
capacity_status=2*ones(n,1); capacity_status(xi<=tol)=0;
capacity_status(xi>=cap-tol)=1;
pa=full_audit.mechanical.phase;
dR=pa.delta_lambda_R(:); dM=pa.delta_lambda_M(:);
composition_mode=get_option(closure,'phase_composition_mode','origin_memory');
sig=struct('key',sprintf('R%s|M%s|C%s|S%s|P%s',mat2str(abs(dR)>tol), ...
    mat2str(abs(dM)>tol),mat2str(capacity_status), ...
    get_option(closure,'phase_solver_mode','lsqnonlin'),composition_mode), ...
    'active_R',abs(dR)>tol,'active_M',abs(dM)>tol, ...
    'capacity_status',capacity_status,'phase_solver_mode', ...
    get_option(closure,'phase_solver_mode','lsqnonlin'), ...
    'phase_composition_mode',composition_mode);
end

function tf=same_signature(a,b)
tf=isstruct(a) && isstruct(b) && isfield(a,'key') && isfield(b,'key') && ...
    strcmp(a.key,b.key);
end

function guess=make_continuation_guess(T,full_audit)
pa=full_audit.mechanical.phase;
plastic_guess=cellfun(@(a)[a.delta_gamma_arp_local; ...
    a.delta_gamma_trip_local],full_audit.mechanical.plastic, ...
    'UniformOutput',false);
guess=struct('temperature_K',T, ...
    'phase',[pa.stress_MPa;pa.delta_lambda_R(:);pa.delta_lambda_M(:)], ...
    'plastic',{plastic_guess});
end

function value=expand_local(value,n)
if isscalar(value), value=repmat(value,n,1); else, value=value(:); end
assert(numel(value)==n && all(isfinite(value)),'Finite per-cell capacity is required.');
end
