function value = hardening_derivative(p, branch, xi, load_memory)
%SMA.HARDENING_DERIVATIVE Branch potential derivative calibrated at endpoints.
%
% This is f_,xi in the Gibbs-force expression. It is constructed so that
% zero-stress forward transformation passes through (Ms,0)->(Mf,1) and
% reverse transformation passes through (As,1)->(Af,0).

xi = min(max(xi,0.0),1.0);
if nargin<4
    load_memory = 0.0;
end
if branch > 0
    surface_temperature = p.Ms - (p.Ms-p.Mf) * xi;
else
    memory = min(max(load_memory,0.0)/p.mt_reverse_memory_scale,1.0)^ ...
        p.mt_reverse_memory_power;
    shift = p.mt_reverse_memory_shift*memory;
    surface_temperature = p.Af+shift - (p.Af-p.As) * xi;
end
value = p.rho_delta_s * (surface_temperature-p.T0) - p.rho_delta_u;
end
