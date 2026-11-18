function [R,ph,pl,tt,d]=candidate_niti_primitive_joint_trial_1d(ctx,value)
% One real source trial; no root solve, state commit or experimental access.
assert(strcmp(ctx.schema_id,'niti_primitive_source_context_v1'));
p=ctx.p;old_phase=ctx.old_phase;old_plastic=ctx.old_plastic;Told=ctx.Told;
sigma_old=ctx.sigma_old;boundary=ctx.boundary;dt=ctx.dt;c=ctx.c;zc=ctx.zc;tc=ctx.tc;
n=ctx.n;m=ctx.m;r0=ctx.r0;w=ctx.w;sR=ctx.sR;sM=ctx.sM;recent=ctx.recent;
capR=ctx.capR;capM=ctx.capM;strain_control=ctx.strain_control;
nr=ctx.nr;nm=ctx.nm;nt=ctx.nt;
value=value(:);assert(numel(value)==numel(ctx.lower) && all(isfinite(value)));
assert(all(value>=ctx.lower) && all(value<=ctx.upper),'Source trial outside declared numerical domain.');
        if strain_control,sigma=value(1);else,sigma=boundary.target;end
        qR=value(nr);qM=value(nm);tt=value(nt);
        dR=sR.*min(dt*r0*qR.^m,capR);dM=sM.*min(dt*r0*qM.^m,capM);
        ph=cell(n,1);pl=cell(n,1);pa=cell(n,1);inc=cell(n,1);
        piR=zeros(n,1);piM=zeros(n,1);deta=zeros(n,1);Y=zeros(n,1);eps=zeros(n,1);
        qresR=zeros(n,1);qresM=zeros(n,1);force=cell(n,1);
        local_z=zc;local_z.force_scalar_hardening=true;
        for j=1:n
            [ph{j},inc{j}]=sma.candidate_zhang_current_mixture_phase_increment_1d(p,old_phase{j},dR(j),dM(j));
            [pl{j},pa{j}]=sma.candidate_arp_trip_inheritance_1d(old_plastic{j},ph{j}.xi_M,sigma,dt,local_z,[]);
            relR=ph{j}.xi_R/max(ph{j}.xi_A+ph{j}.xi_R,realmin);
            piR(j)=p.R_H*sigma+p.R_entropy*(tt(j)-c.Tref_K)-p.R_tau0-p.R_hardening*relR;
            NB=sign(old_plastic{j}.eps_p_AR);if NB==0,NB=sM(j);end
            Btr=pl{j}.Omega*NB;
            source=sma.candidate_zhang_source_force_representation_1d( ...
                inc{j}.g_mixed,sigma,Btr,pl{j}.Omega,c.source_force_representation);
            deta(j)=sma.mixed_mt_entropy(p,inc{j}.weight_R,pl{j}.DeltaEtaP_MPaK);
            tau=sma.candidate_zhang_source_mt_resistance(p,sM(j),ph{j}.xi_M,recent(j),pl{j}.xi_Mr);
            [Y(j),ya]=sma.candidate_zhang_kinetic_threshold_1d(p,recent(j),pl{j}.Yp);
            elastic=.5*sigma^2*(1/p.EM-1/p.EA);thermal=deta(j)*(tt(j)-c.Tref_K);
            piM(j)=source.mechanical_force_MPa+elastic+thermal-tau;
            availR=(sR(j)>0 && old_phase{j}.xi_A>0) || (sR(j)<0 && old_phase{j}.xi_R>0);
            availM=(sM(j)>0 && capM(j)>0) || (sM(j)<0 && old_phase{j}.xi_M>old_plastic{j}.xi_Mr);
            qresR(j)=sma.candidate_source_q_kkt_residual_1d(piR(j),c.Y_R_MPa,sR(j),qR(j),dt,r0,m,capR(j),availR);
            qresM(j)=sma.candidate_source_q_kkt_residual_1d(piM(j),Y(j),sM(j),qM(j),dt,r0,m,capM(j),availM);
            P=ph{j}.xi_M*pl{j}.eps_p_M+(1-ph{j}.xi_M)*pl{j}.eps_p_AR;
            eps(j)=sigma/sma.effective_modulus(p,ph{j}.xi_M)+p.alphaA*(tt(j)-p.T0)+ph{j}.eps_RT+ph{j}.eps_MT+P;
            force{j}=struct('mechanical_MPa',source.mechanical_force_MPa,'elastic_MPa',elastic, ...
                'thermal_MPa',thermal,'tau_MPa',tau,'pi_MPa',piM(j),'Y_MPa',Y(j), ...
                'Mr_used',pl{j}.xi_Mr,'Omega_used',pl{j}.Omega,'Yp_used',pl{j}.Yp, ...
                'DeltaEtaP_used',pl{j}.DeltaEtaP_MPaK,'recent',recent(j),'threshold',ya);
        end
        phase_audit=struct('stress_MPa',sigma,'delta_lambda_R',dR,'delta_lambda_M',dM, ...
            'terms',struct('pi_R_MPa',piR,'pi_M_MPa',piM,'delta_eta_mix_MPaK',deta));
        [Tbalance,heat]=sma.candidate_gradient_heat_balance_1d(p,Told,tt,sigma_old,phase_audit,pa,dt,tc);
        R=[qresR;qresM;Tbalance-tt];macro=w.'*eps;
        if strain_control,R=[(macro-boundary.target)/p.Hmax_i;R];end
        d=struct('stress_MPa',sigma,'macro_strain',macro,'phase',phase_audit,'plastic_audit',{pa}, ...
            'force_ledger',{force},'heat',heat,'q_R',qR,'q_M',qM,'q_KKT_R',qresR,'q_KKT_M',qresM, ...
            'source_same_trial_history_from_accepted_candidate',true,'local_plastic_elimination','exact_shared_hardening_scalar');
end
