function [phase,plastic,T,audit]=candidate_niti_primitive_joint_step_1d( ...
    p,old_phase,old_plastic,Told,sigma_old,boundary,dt,c,zc,tc,guess)
% Source-only ell=0 primitive joint residual. No outer history Picard.
% Unknowns: stress for strain control, both natural-q coordinates, local T.
% Local plasticity and source histories are exact candidate eliminations.
ctx=sma.candidate_niti_primitive_joint_context_1d(p,old_phase,old_plastic,Told,sigma_old,boundary,dt,c,zc,tc);
n=ctx.n;m=ctx.m;r0=ctx.r0;strain_control=ctx.strain_control;
offset=ctx.offset;nr=ctx.nr;nm=ctx.nm;nt=ctx.nt;lower=ctx.lower;upper=ctx.upper;
qcapR=ctx.qcapR;qcapM=ctx.qcapM;
v=zeros(offset+3*n,1);v(nt)=Told;
if strain_control,v(1)=sigma_old;end
if nargin>=11 && ~isempty(guess)
    if strain_control && isfield(guess,'sigma_MPa'),v(1)=guess.sigma_MPa;end
    if isfield(guess,'temperature_K'),v(nt)=guess.temperature_K(:);end
    if isfield(guess,'delta_lambda_R'),v(nr)=(abs(guess.delta_lambda_R(:))/(dt*r0)).^(1/m);end
    if isfield(guess,'delta_lambda_M'),v(nm)=(abs(guess.delta_lambda_M(:))/(dt*r0)).^(1/m);end
end
v=min(max(v,lower),upper);active=find(upper>lower);prototype=v;
calls=0;last_failure='';
options=optimoptions('lsqnonlin','Display','off','FiniteDifferenceType','central', ...
    'ScaleProblem','jacobian','FunctionTolerance',1e-12,'StepTolerance',1e-12, ...
    'OptimalityTolerance',1e-10,'MaxIterations',100,'MaxFunctionEvaluations',2500);
[x,~,~,flag,output]=lsqnonlin(@reduced,v(active),lower(active),upper(active),options);
v(active)=x;[R,phase,plastic,T,detail]=evaluate(v);
% Reflective box least squares only approaches an active q bound from
% inside. Its optimizer tolerance does not enforce the source q gate.
% A projected Newton correction solves the SAME coupled residual and may
% land exactly on that inventory face. Every trial re-evaluates phase,
% plastic histories, thermal balance and the imposed control together.
newton_trace=[];
for correction=1:15
    if norm(R(active),inf)<=1e-9,break;end
    J=zeros(numel(active));
    for col=1:numel(active)
        index=active(col);h=2e-6*max(1,abs(v(index)));
        vp=v;vm=v;vp(index)=min(v(index)+h,upper(index));
        vm(index)=max(v(index)-h,lower(index));
        Rp=evaluate(vp);Rm=evaluate(vm);
        J(:,col)=(Rp(active)-Rm(active))/(vp(index)-vm(index));
    end
    if rcond(J)<1e-14,break;end
    increment=-J\R(active);before=norm(R(active),2);accepted=false;
    for line=0:16
        trial=v;trial(active)=min(max(v(active)+(2^-line)*increment,lower(active)),upper(active));
        try
            Rt=evaluate(trial);
            if norm(Rt(active),2)<before || norm(Rt(active),inf)<=1e-9
                v=trial;R=Rt;accepted=true;break;
            end
        catch err
            last_failure=[err.identifier ': ' err.message];
        end
    end
    newton_trace=[newton_trace;before,norm(R(active),2),double(accepted)]; %#ok<AGROW>
    if ~accepted,break;end
end
[R,phase,plastic,T,detail]=evaluate(v);
qerror=max(abs(R(offset+(1:2*n))));terror=max(abs(R(nt)));
if strain_control,control_error=abs(detail.macro_strain-boundary.target);else,control_error=0;end
assert(qerror<=2e-8 && terror<=1e-7 && control_error<=2e-9, ...
    ['Primitive source residual failed: q %.9g T %.9g strain %.9g flag %d; ' ...
    'qR=%s qM=%s rR=%s rM=%s capR=%s capM=%s; %s'], ...
    qerror,terror,control_error,flag,mat2str(v(nr),15),mat2str(v(nm),15), ...
    mat2str(R(nr),15),mat2str(R(nm),15),mat2str(qcapR,15),mat2str(qcapM,15),last_failure);
audit=detail;audit.q_residual=qerror;audit.thermal_residual_K=terror;
audit.strain_control_residual=detail.macro_strain-boundary.target;
if ~strain_control,audit.strain_control_residual=NaN;end
audit.iterations=output.iterations;audit.function_calls=calls;audit.exitflag=flag;
audit.least_squares_exitflag=flag;audit.projected_newton_steps=size(newton_trace,1);
audit.projected_newton_residual_trace=newton_trace;
audit.solver_owner='primitive_source_q_phase_temperature_with_local_history_elimination';
audit.effective_identity=struct('strict_source_kinetics',true, ...
    'candidate_history_owner','same_trial_full_history_from_candidate_M_and_stress', ...
    'recent_history_owner','committed_cumulative_max_within_trial', ...
    'thermal_owner','per_node_source_heat_balance_no_homogeneous_reduction', ...
    'phase_solver_mode',c.phase_solver_mode,'local_plastic_solver','exact_shared_hardening_scalar');
audit.outer_history_fixed_point_used=false;audit.material_qualified=false;
audit.production_connected=false;audit.observations_accessed=false;audit.experimental_backfill=false;

    function R=reduced(x)
        trial=prototype;trial(active)=x;
        try
            R=evaluate(trial);
        catch err
            last_failure=[err.identifier ': ' err.message];
            R=1e4*ones(numel(prototype),1);
        end
    end
    function [R,ph,pl,tt,d]=evaluate(value)
        calls=calls+1;
        [R,ph,pl,tt,d]=sma.candidate_niti_primitive_joint_trial_1d(ctx,value);
    end
end
