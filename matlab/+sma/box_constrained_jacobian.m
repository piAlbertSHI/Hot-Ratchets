function J = box_constrained_jacobian(fun, x, relative_step, lower, upper)
%SMA.BOX_CONSTRAINED_JACOBIAN Finite differences inside feasible bounds.

x = x(:);
lower = lower(:);
upper = upper(:);
assert(isequal(size(x),size(lower),size(upper)), ...
    'State and bound vectors must have identical sizes.');
assert(all(lower<=x) && all(x<=upper), ...
    'The Jacobian evaluation point must lie inside its bounds.');

f0 = fun(x);
J = zeros(numel(f0),numel(x));
for k = 1:numel(x)
    step = relative_step*max(1.0,abs(x(k)));
    forward_room = upper(k)-x(k);
    backward_room = x(k)-lower(k);
    if backward_room < step && forward_room > 0.0
        h = min(step,forward_room);
        xp = x;
        xp(k) = xp(k)+h;
        J(:,k) = (fun(xp)-f0)/h;
    elseif forward_room < step && backward_room > 0.0
        h = min(step,backward_room);
        xm = x;
        xm(k) = xm(k)-h;
        J(:,k) = (f0-fun(xm))/h;
    else
        h = min([step,forward_room,backward_room]);
        assert(h>0.0, ...
            'A fixed variable cannot be differentiated numerically.');
        xp = x;
        xm = x;
        xp(k) = xp(k)+h;
        xm(k) = xm(k)-h;
        J(:,k) = (fun(xp)-fun(xm))/(2.0*h);
    end
end
end
