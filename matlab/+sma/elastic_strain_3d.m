function h_e = elastic_strain_3d(p, xi, tau)
%SMA.ELASTIC_STRAIN_3D Isotropic compliance mixture.

E = sma.effective_modulus(p, xi);
nu = sma.effective_poisson(p, xi);
h_e = (1.0 + nu) / E * tau - nu / E * trace(tau) * eye(3);
end
