function [pi_value,phi,terms] = driving_force_1d( ...
    p,branch,sigma,beta,T,xi,lambda_tr,lambda_tp,zeta_d,xi_recent, ...
    training_activity,load_memory,population_activity,direct_Yp,direct_Omega, ...
    delta_eta_AR_M)
%SMA.DRIVING_FORCE_1D Full reduced thermodynamic driving force.

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
delta_compliance = 1.0/p.EM - 1.0/p.EA;
delta_alpha = p.alphaM-p.alphaA;
hardening = sma.hardening_derivative(p,branch,xi,load_memory);
thermal_capacity = -p.rho_delta_c * ...
    ((T-p.T0)-T*log(max(T,eps)/p.T0));
terms = struct();
transformation_direction=sign(lambda_tr);
if transformation_direction==0, transformation_direction=1.0; end
direct_internal_stress=p.zhang_direct_Omega_coupling* ...
    max(direct_Omega,0.0)*transformation_direction;
terms.transformation_work = ...
    (sigma+beta+direct_internal_stress)*lambda_tr;
terms.direct_internal_stress=direct_internal_stress;
terms.trip_work = sigma*lambda_tp;
terms.compliance = 0.5*sigma^2*delta_compliance;
terms.thermal_expansion = sigma*delta_alpha*(T-p.T0);
terms.entropy = delta_eta_AR_M*T;
terms.heat_capacity = thermal_capacity;
terms.internal_energy = -p.rho_delta_u;
terms.hardening = -hardening;
% The T0 reference is subtracted from the entropy term so the fitted
% internal-energy/hardening constants remain numerically well scaled.
terms.entropy_difference = delta_eta_AR_M;
terms.entropy = delta_eta_AR_M*(T-p.T0);
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
