function alpha = thermal_expansion(p, xi)
%SMA.THERMAL_EXPANSION Linear phase mixture.

xi = min(max(xi, 0.0), 1.0);
alpha = p.alphaA + xi * (p.alphaM - p.alphaA);
end
