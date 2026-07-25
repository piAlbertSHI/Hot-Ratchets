function [Hcur, Hmax] = hcur(p, equivalent_stress, zeta_d)
%SMA.HCUR Current and cycle-evolved transformation strain magnitude.

degradation_activity = max(zeta_d-p.zeta_H0,0.0);
Hmax = p.Hmax_f + (p.Hmax_i - p.Hmax_f) * ...
    exp(-p.lambda_H * degradation_activity);
activation = max(equivalent_stress - p.sigma_crit, 0.0);
Hcur = p.Hmin + (Hmax - p.Hmin) * (1.0 - exp(-p.kt * activation));
Hcur = min(max(Hcur, p.Hmin), Hmax);
end
