function value = transformation_resistance( ...
    p,branch,zeta_d,xi,xi_recent,training_activity,load_memory, ...
    population_activity,direct_Yp)
%SMA.TRANSFORMATION_RESISTANCE Cycle-evolved forward/reverse resistance.

if nargin<4
    xi = 0.0;
end
if nargin<5
    xi_recent = max(xi,p.numerics.xi_tol);
end
if nargin<6
    training_activity = zeta_d;
end
if nargin<7
    load_memory = 0.0;
end
if nargin<8 || isempty(population_activity)
    population_activity=zeros(p.Y_activity_population_bins,1);
end
if nargin<9 || isempty(direct_Yp)
    direct_Yp=0.0;
end
assert(~(p.enable_zhang_source_mt_resistance && ...
    p.zhang_direct_Yp_coupling>0), ...
    'sma:transformation_resistance:ConflictingMTOwners', ...
    ['Zhang source MT resistance and direct Yp coupling are mutually ' ...
    'exclusive; enabling both would count the MT surface twice.']);
% The recent transformed fraction is a shared state argument for both the
% dissipation-evolution and nonlinear-hardening closures.  Define it before
% either optional mechanism so the mechanisms remain independently
% switchable in the unified model.
xi_recent_eff = min(max([xi,xi_recent,0.0]),1.0);
activity_source = zeta_d;
if p.Y_activity_use_population_path
    % zeta_p is committed only after a completed reverse traversal.  It is
    % therefore a cycle-frozen first moment of the active transformation
    % population, not a cycle-number lookup or an observation override.
    activity_source = training_activity;
end
activity = max(activity_source-p.zeta_Y0,0.0);
legacy_factor = sma.activity_relaxation_factor(p,activity,xi_recent_eff);
population_factor=sma.population_relaxation_factor( ...
    p,population_activity,xi);
mix=p.Y_activity_population_mix;
factor=(1.0-mix)*legacy_factor+mix*population_factor;
high_activity = min(max(zeta_d,0.0)/p.mt_Yp_activity_scale,1.0)^ ...
    p.mt_Yp_activity_power;
memory_activity = min(max(load_memory,0.0)/p.mt_Yp_memory_scale,1.0)^ ...
    p.mt_Yp_memory_power;
if branch>0
    value = p.Yf_f+(p.Yf-p.Yf_f)*factor;
else
    value = p.Yr_f+(p.Yr-p.Yr_f)*factor;
end
if p.enable_zhang_source_mt_resistance
    value=sma.candidate_zhang_source_mt_resistance( ...
        p,branch,xi,xi_recent_eff,load_memory);
end
if p.enable_mt_dissipation_evolution && ~p.enable_zhang_source_mt_resistance
    Yp_activity = max(training_activity,0.0);
    if p.mt_Yp_use_cumulative_activity
        Yp_activity = max(zeta_d,0.0);
    end
    Yp_saturation = p.mt_Yp_sat+ ...
        p.mt_Yp_high_sat*memory_activity;
    Yp = Yp_saturation*(1.0-exp( ...
        -p.mt_Yp_rate*Yp_activity));
    if branch>0
        Yf_amplitude = p.mt_Y_d*xi_recent_eff^p.mt_Y_n;
        value = value+Yf_amplitude-Yp;
    else
        % Yr is only the explicit correction left after the As/Af endpoint
        % hardening has been closed. Yp therefore acts on that total
        % reduced surface and is not clipped at Yr=0.
        reverse_memory = min(max(load_memory,0.0)/ ...
            p.mt_Yp_reverse_memory_scale,1.0)^ ...
            p.mt_Yp_reverse_memory_power;
        reverse_training_factor = exp(-p.mt_Yp_reverse_training_decay* ...
            max(training_activity,0.0));
        value = value-p.mt_Yp_reverse_scale*reverse_training_factor* ...
            (1.0-reverse_memory)*Yp;
    end
end
if p.enable_mt_nonlinear_hardening && ~p.enable_zhang_source_mt_resistance
    xi_eff = min(max(xi,0.0),1.0);
    if branch>0
        training_factor = min(max(training_activity,0.0)/ ...
            p.mt_hfor_training_scale,1.0)^p.mt_hfor_training_power;
        transition_progress = min( ...
            xi_eff/p.mt_hfor_transition_width,1.0);
        transition_dip = p.mt_hfor_transition_dip* ...
            4.0*transition_progress*(1.0-transition_progress);
        shoulder_progress = min(max( ...
            (xi_eff-p.mt_hfor_shoulder_start)/ ...
            p.mt_hfor_shoulder_width,0.0),1.0);
        shoulder_memory = min(max(load_memory,0.0)/ ...
            p.mt_hfor_shoulder_memory_scale,1.0)^ ...
            p.mt_hfor_shoulder_memory_power;
        shoulder = p.mt_hfor_shoulder* ...
            (1.0-shoulder_memory)* ...
            shoulder_progress^2*(3.0-2.0*shoulder_progress)* ...
            min(max(training_activity,0.0)/ ...
            p.mt_hfor_shoulder_training_scale,1.0)^ ...
            p.mt_hfor_shoulder_training_power;
        hfor = p.mt_hfor+p.mt_hfor_training_extra*training_factor;
        nonlinear = hfor*xi_eff^p.mt_n+ ...
            p.mt_hfor_activity*high_activity* ...
            xi_eff^p.mt_hfor_activity_power+shoulder-transition_dip;
    else
        reverse_progress = 0.0;
        if xi_recent>p.numerics.xi_tol
            reverse_progress = 1.0-xi_eff/xi_recent;
        end
        reverse_progress = min(max(reverse_progress,0.0),1.0);
        reverse_span = p.mt_hrev_span* ...
            (1.0-xi_recent_eff)^p.mt_hrev_span_load_power* ...
            (2.0*reverse_progress-1.0);
        reverse_finish_low = p.mt_hrev_finish_low* ...
            (1.0-xi_recent_eff)^p.mt_hrev_finish_low_load_power* ...
            reverse_progress^p.mt_hrev_finish_low_power* ...
            max(0.0,1.0-min(max(training_activity,0.0)/ ...
            p.mt_hrev_finish_low_training_scale,1.0)^ ...
            p.mt_hrev_finish_low_training_power);
        nonlinear = p.mt_hrev*xi_eff^p.mt_n+ ...
            p.mt_hrev_finish*xi_recent_eff^ ...
            p.mt_hrev_finish_load_power* ...
            reverse_progress^p.mt_hrev_finish_power+reverse_span+ ...
            reverse_finish_low;
    end
    value = value+p.mt_H*xi_eff+nonlinear;
end
if p.zhang_direct_Yp_coupling>0
    % Zhang Eqs. (54)-(58): Y = Yf0 + d*(xi_recent)^n - Yp.
    % The direct candidate is a replacement for the reduced mt_Yp closure;
    % validation therefore disables enable_mt_dissipation_evolution to
    % avoid counting d or Yp twice.
    direct_value=p.mt_Y_d*xi_recent_eff^p.mt_Y_n-max(direct_Yp,0.0);
    value=value+p.zhang_direct_Yp_coupling*direct_value;
end
end
