function ctx=candidate_niti_primitive_joint_context_1d(p,old_phase,old_plastic,Told,sigma_old,boundary,dt,c,zc,tc)
% Exact source primitive scope/bounds, shared by solver and derivative probes.
assert(isfield(c,'strict_source_kinetics') && isequal(c.strict_source_kinetics,true), ...
    'sma:primitive_source:KineticIdentity','Primitive source solver requires explicit strict-q identity.');
assert(isfield(c,'residual_martensite_coupling') && ...
    strcmp(c.residual_martensite_coupling,'same_trial_full_history') && ...
    isfield(c,'recent_M_owner') && strcmp(c.recent_M_owner,'cumulative_max') && ...
    isfield(c,'same_trial_history_owner') && strcmp(c.same_trial_history_owner,'committed_source'), ...
    'sma:primitive_source:HistoryIdentity', ...
    'Candidate histories must be full same-trial; recent M remains committed cumulative maximum within trial.');
assert(isfield(c,'homogeneous_thermal_reduction') && isequal(c.homogeneous_thermal_reduction,false), ...
    'sma:primitive_source:ThermalIdentity','Primitive source solver retains every nodal temperature.');
assert(isfield(c,'phase_solver_mode') && strcmp(c.phase_solver_mode,'primitive_joint_projected_newton'), ...
    'sma:primitive_source:SolverIdentity','Explicit primitive joint solver identity is required.');
assert(~isfield(c,'population_count') && ...
    (~isfield(c,'resistance_multiplier') || all(c.resistance_multiplier(:)==1)), ...
    'sma:primitive_source:UnsupportedClosure','Population kinetics or resistance scaling is outside this primitive scope.');
for name={'transformation_exponent','transformation_exponent_R','transformation_exponent_M'}
    if isfield(c,name{1}),assert(c.(name{1})==p.zhang_m_tr, ...
        'sma:primitive_source:KineticIdentity','Kinetic exponent override differs from the bound source card.');end
end
if isfield(c,'transformation_reference_rate_per_s')
    assert(c.transformation_reference_rate_per_s==p.zhang_reference_rate_per_s, ...
        'sma:primitive_source:KineticIdentity','Reference rate override differs from the bound source card.');
end
assert(p.case_id==1 && c.kappa_MPa_length2==0 && dt>0 && ...
    strcmp(c.phase_composition_mode,'source_current_mixture') && ...
    strcmp(zc.residual_m_event_policy,'event_boundary'), ...
    'Primitive kernel currently binds NiTi source ell=0/current-mixture/untruncated-event only.');
assert(~isfield(c,'enable_path_accommodation') || ~c.enable_path_accommodation);
assert(~isfield(c,'enable_thermal_accommodation') || ~c.enable_thermal_accommodation);
assert(all(c.resistance_shift_MPa(:)==0) && all(c.f_transformable(:)==1), ...
    'Noncanonical resistance or transformable capacity needs its own residual owner.');
assert(~isfield(c,'eq39_ablation') || ~c.eq39_ablation.enabled);
n=numel(old_phase);assert(n>=2 && numel(old_plastic)==n);
Told=Told(:);w=c.weights(:);assert(numel(w)==n && all(w>0) && abs(sum(w)-1)<1e-12);
sR=expand(c.direction_R,n);sM=expand(c.direction_M,n);recent=expand(c.xi_recent,n);
assert(all(abs(sR)==1) && all(abs(sM)==1) && all(Told>40));
m=p.zhang_m_tr;r0=p.zhang_reference_rate_per_s;
capR=zeros(n,1);capM=zeros(n,1);
for i=1:n
    assert(old_phase{i}.xi_A+old_phase{i}.xi_R>0,'Pure-M parent needs a separate source composition chart.');
    if sR(i)>0,capR(i)=old_phase{i}.xi_A;else,capR(i)=old_phase{i}.xi_R;end
    if sM(i)>0,capM(i)=old_phase{i}.xi_A+old_phase{i}.xi_R;
    else,capM(i)=max(old_phase{i}.xi_M-event_lower(old_phase{i},old_plastic{i},zc),0);end
end
qcapR=(capR/(dt*r0)).^(1/m);qcapM=(capM/(dt*r0)).^(1/m);
strain_control=strcmp(boundary.kind,'strain');assert(strain_control || strcmp(boundary.kind,'stress'));
offset=double(strain_control);nr=offset+(1:n);nm=offset+n+(1:n);nt=offset+2*n+(1:n);
lower=[zeros(offset+2*n,1);Told-40];upper=[zeros(offset,1);qcapR;qcapM;Told+40];
if strain_control,lower(1)=-p.numerics.stress_solve_bound;upper(1)=p.numerics.stress_solve_bound;end
ctx=struct('schema_id','niti_primitive_source_context_v1','p',p, ...
    'old_phase',{old_phase},'old_plastic',{old_plastic},'Told',Told, ...
    'sigma_old',sigma_old,'boundary',boundary,'dt',dt,'c',c,'zc',zc,'tc',tc, ...
    'n',n,'m',m,'r0',r0,'w',w,'sR',sR,'sM',sM,'recent',recent, ...
    'capR',capR,'capM',capM,'qcapR',qcapR,'qcapM',qcapM, ...
    'strain_control',strain_control,'offset',offset,'nr',nr,'nm',nm,'nt',nt, ...
    'lower',lower,'upper',upper);
end

function value=expand(value,n)
if isscalar(value),value=repmat(value,n,1);else,value=value(:);end
assert(numel(value)==n && all(isfinite(value)));
end
function r=event_lower(ph,pl,z)
r=pl.xi_Mr;
if ph.xi_M<=r+2e-12,return;end
F=z.zhang_cxi_sqrt_m*sqrt(max(pl.rho_M,0));a=z.zhang_a1;
fun=@(M)M-F-(pl.xi_Mr-F)*exp(-a*(ph.xi_M-M));
if fun(0)>=0 || fun(ph.xi_M)<=0,return;end
lo=0;hi=ph.xi_M;
for k=1:90
    mid=(lo+hi)/2;if fun(mid)>0,hi=mid;else,lo=mid;end
end
r=(lo+hi)/2;
end
