function J = numerical_jacobian(fun, x, relative_step)
%SMA.NUMERICAL_JACOBIAN Central-difference Jacobian for local residuals.

f0 = fun(x);
J = zeros(numel(f0), numel(x));
for k = 1:numel(x)
    step = relative_step * max(1.0, abs(x(k)));
    xp = x;
    xm = x;
    xp(k) = xp(k) + step;
    xm(k) = xm(k) - step;
    J(:, k) = (fun(xp) - fun(xm)) / (2.0 * step);
end
end
