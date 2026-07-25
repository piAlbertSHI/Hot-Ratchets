function E = effective_modulus(p, xi)
%SMA.EFFECTIVE_MODULUS Reuss mixture consistent with S(xi).

xi = min(max(xi, 0.0), 1.0);
E = 1.0 / ((1.0 - xi) / p.EA + xi / p.EM);
end
