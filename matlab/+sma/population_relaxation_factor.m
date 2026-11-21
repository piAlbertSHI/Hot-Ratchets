function factor = population_relaxation_factor(p,activity,xi)
%SMA.POPULATION_RELAXATION_FACTOR Marginal threshold-population resistance.

activity=activity(:);
n=p.Y_activity_population_bins;
if numel(activity)~=n || any(~isfinite(activity)) || any(activity<0)
    error('sma:population_relaxation_factor:InvalidState', ...
        'Population activity is invalid.');
end
q=((1:n).'-0.5)/n;
local=p.Y_activity_fast_weight*exp(-p.Y_activity_fast_rate*activity)+ ...
    (1.0-p.Y_activity_fast_weight)*exp(-p.Y_activity_slow_rate*activity);
kernel=exp(-0.5*((q-min(max(xi,0.0),1.0))/ ...
    p.Y_activity_population_width).^2);
weights=kernel/sum(kernel);
factor=sum(weights.*local);
factor=min(max(factor,0.0),1.0);
end
