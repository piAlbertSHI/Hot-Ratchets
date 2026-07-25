function A = voigt_to_tensor(v)
%SMA.VOIGT_TO_TENSOR Tensorial-shear Voigt ordering 11,22,33,12,23,13.

v = v(:);
assert(numel(v) == 6, 'Voigt vector must contain six entries.');
A = [v(1), v(4), v(6); v(4), v(2), v(5); v(6), v(5), v(3)];
end
