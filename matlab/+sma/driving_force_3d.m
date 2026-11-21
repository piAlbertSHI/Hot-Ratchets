function [pi_value,phi,terms] = driving_force_3d( ...
    p,branch,tau,beta,T,xi,Lambda,Lambda_tp,zeta_d,xi_recent, ...
    training_activity,load_memory,population_activity,direct_Yp,direct_Omega, ...
    delta_eta_AR_M)
%SMA.DRIVING_FORCE_3D Full tensor thermodynamic driving force.

if nargin<9
    zeta_d = 0.0;
end
if nargin<10
    xi_recent = max(xi,p.numerics.xi_tol);
end
if nargin<11
    training_activity = zeta_d;
end
if nargin<12
    load_memory = 0.0;
end
if nargin<13
    population_activity=zeros(p.Y_activity_population_bins,1);
end
if nargin<14, direct_Yp=0.0; end
if nargin<15, direct_Omega=0.0; end
if nargin<16, delta_eta_AR_M=p.rho_delta_s; end
delta_strain = sma.phase_elastic_strain_3d(p.EM,p.nuM,tau) - ...
    sma.phase_elastic_strain_3d(p.EA,p.nuA,tau);
delta_alpha = p.alphaM-p.alphaA;
hardening = sma.hardening_derivative(p,branch,xi,load_memory);
thermal_capacity = -p.rho_delta_c * ...
    ((T-p.T0)-T*log(max(T,eps)/p.T0));
terms = struct();
Lambda_eq=sqrt(max((2.0/3.0)*sum(Lambda.^2,'all'),0.0));
if Lambda_eq>eps
    Ntr=Lambda/Lambda_eq;
else
    Ntr=zeros(3);
end
direct_internal_stress=(2.0/3.0)*p.zhang_direct_Omega_coupling* ...
    max(direct_Omega,0.0)*Ntr;
terms.transformation_work = ...
    sum((tau+beta+direct_internal_stress).*Lambda,'all');
terms.direct_internal_stress=direct_internal_stress;
terms.trip_work = sum(tau.*Lambda_tp,'all');
terms.compliance = 0.5*sum(tau.*delta_strain,'all');
terms.thermal_expansion = trace(tau)*delta_alpha*(T-p.T0);
terms.entropy_difference = delta_eta_AR_M;
terms.entropy = delta_eta_AR_M*(T-p.T0);
terms.heat_capacity = thermal_capacity;
terms.internal_energy = -p.rho_delta_u;
terms.hardening = -hardening;
pi_value = terms.transformation_work+terms.trip_work+terms.compliance+ ...
    terms.thermal_expansion+terms.entropy+terms.heat_capacity+ ...
    terms.internal_energy+terms.hardening;
resistance = sma.transformation_resistance( ...
    p,branch,zeta_d,xi,xi_recent,training_activity,load_memory, ...
    population_activity,direct_Yp);
terms.resistance = resistance;
if branch > 0
    phi = pi_value-resistance;
else
    phi = -pi_value-resistance;
end
end
