function v = tensor_to_voigt(A)
%SMA.TENSOR_TO_VOIGT Tensorial-shear Voigt ordering 11,22,33,12,23,13.

v = [A(1,1); A(2,2); A(3,3); A(1,2); A(2,3); A(1,3)];
end
