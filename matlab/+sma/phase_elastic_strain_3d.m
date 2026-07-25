function h_e = phase_elastic_strain_3d(E,nu,tau)
%SMA.PHASE_ELASTIC_STRAIN_3D Isotropic phase compliance action.

h_e = (1.0+nu)/E*tau - nu/E*trace(tau)*eye(3);
end
