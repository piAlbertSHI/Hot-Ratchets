function [next,audit] = zhang_history_evolution_direct(previous,xi_new, ...
    delta_eps_arp,delta_eps_tp)
%SMA.ZHANG_HISTORY_EVOLUTION_DIRECT Zhang 2022 Eqs. (48)-(60) shadow map.
% All coefficients are the direct Table 2 values.  This evaluator is
% deliberately isolated from the production return map and never reads an
% observation or a cycle number.

assert(isstruct(previous),'Previous Zhang history must be a struct.');
required={'xi','xi_Mr','rho_M','rho_AR','D_M','D_AR','DeltaEtaP', ...
    'Omega','Yp','Lambda_M'};
for k=1:numel(required)
    assert(isfield(previous,required{k}), ...
        'sma:ZhangDirect:MissingState','Missing field %s.',required{k});
end
assert(isfinite(xi_new) && xi_new>=0 && xi_new<=1, ...
    'sma:ZhangDirect:PhaseBounds','xi_new must be in [0,1].');
assert(isfinite(delta_eps_arp) && isfinite(delta_eps_tp), ...
    'sma:ZhangDirect:Increment','Plastic increments must be finite.');
validate_state(previous);

% Zhang Table 2 direct values.
k1=2.0e9; k2=50.0; a1=0.6; a2=0.1;
cOmega=92.0; cxi=3.5e-8; etaM=1.23e-9; etaA=1.0e-9; cY=4.5;
phase_tol=1.0e-12;
xi_old=previous.xi;
dxi=xi_new-xi_old;
dxi_plus=max(dxi,0.0); dxi_minus=max(-dxi,0.0);
dLambda=abs(dxi);

% Eq. (48-b), integrated through q=sqrt(rho_AR), including the
% non-Lipschitz zero-density start without inserting a fitted seed.
gamma_increment=abs(delta_eps_arp)+abs(delta_eps_tp);
q_old=sqrt(max(previous.rho_AR,0.0));
if k2>0
    q_sat=k1/k2;
    q_new=q_sat+(q_old-q_sat)*exp(-0.5*k2*gamma_increment);
else
    q_new=q_old+0.5*k1*gamma_increment;
end
rho_ar_star=max(q_new,0.0)^2;

% Eqs. (48-a), (49-b), and forward/reverse inheritance.
if dxi_plus>0
    D_M=previous.D_M+dxi_plus*rho_ar_star;
    D_AR=max(1.0-xi_new,0.0)*rho_ar_star;
elseif dxi_minus>0
    D_M=max(xi_new,0.0)*max(previous.rho_M,0.0);
    D_AR=max(1.0-xi_new,0.0)*rho_ar_star+ ...
        dxi_minus*max(previous.rho_M,0.0);
else
    D_M=previous.D_M; D_AR=previous.D_AR;
end
D_M=max(D_M,0.0); D_AR=max(D_AR,0.0);
if xi_new>phase_tol
    rho_M=D_M/xi_new;
else
    rho_M=max(previous.rho_M,0.0);
end
if 1.0-xi_new>phase_tol
    rho_AR=D_AR/(1.0-xi_new);
else
    rho_AR=max(rho_ar_star,previous.rho_AR);
end

% Eqs. (53-a), (53-b), (55)-(58): c_Omega, c_xi, and c_Y are direct
% constants per Eq. (57), while residual martensite depends on rho_M.
xi_target=min(max(cxi*sqrt(rho_M),0.0),xi_new);
xi_Mr=xi_target+(previous.xi_Mr-xi_target)*exp(-a1*dLambda);
xi_Mr=min(max(xi_Mr,0.0),xi_new);
Omega=cOmega+(previous.Omega-cOmega)*exp(-a2*dLambda);
Yp=cY+(previous.Yp-cY)*exp(-a2*dLambda);

% Eq. (60), dislocation-induced configurational entropy increment.
DeltaEtaP=previous.DeltaEtaP+etaM*(rho_M-previous.rho_M)- ...
    etaA*(rho_AR-previous.rho_AR);
next=struct('xi',xi_new,'xi_Mr',xi_Mr,'rho_M',rho_M, ...
    'rho_AR',rho_AR,'D_M',D_M,'D_AR',D_AR,'DeltaEtaP',DeltaEtaP, ...
    'Omega',Omega,'Yp',Yp,'Lambda_M',previous.Lambda_M+dLambda);
if nargout>1
    rho_mix=xi_new*rho_M+(1.0-xi_new)*rho_AR;
    audit=struct('dxi',dxi,'dLambda_M',dLambda, ...
        'gamma_increment',gamma_increment,'rho_mix',rho_mix, ...
        'D_balance_error',D_M+D_AR-rho_mix, ...
        'xi_inventory',xi_new-xi_Mr,'state_finite', ...
        all(isfinite([xi_Mr,rho_M,rho_AR,D_M,D_AR,DeltaEtaP, ...
        Omega,Yp,next.Lambda_M])));
end
end

function validate_state(s)
v=[s.xi,s.xi_Mr,s.rho_M,s.rho_AR,s.D_M,s.D_AR,s.DeltaEtaP, ...
    s.Omega,s.Yp,s.Lambda_M];
assert(all(isfinite(v)),'sma:ZhangDirect:NonfiniteState', ...
    'Zhang shadow state is nonfinite.');
assert(s.xi>=0 && s.xi<=1 && s.xi_Mr>=0 && ...
    s.xi_Mr<=s.xi+1e-12 && s.rho_M>=0 && s.rho_AR>=0 && ...
    s.D_M>=0 && s.D_AR>=0 && s.Lambda_M>=0, ...
    'sma:ZhangDirect:Bounds','Zhang shadow state bounds invalid.');
end
