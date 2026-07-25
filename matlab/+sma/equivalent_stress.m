function [seq, dev] = equivalent_stress(tau)
%SMA.EQUIVALENT_STRESS J2 equivalent stress and deviator.

dev = tau - trace(tau) / 3.0 * eye(3);
seq = sqrt(max(1.5 * sum(dev(:) .* dev(:)), 0.0));
end
