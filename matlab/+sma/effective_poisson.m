function nu = effective_poisson(p, xi)
%SMA.EFFECTIVE_POISSON Linear phase mixture for Poisson ratio.

xi = min(max(xi, 0.0), 1.0);
nu = p.nuA + xi * (p.nuM - p.nuA);
end
