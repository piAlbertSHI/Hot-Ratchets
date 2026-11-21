function activity = update_population_activity(p,previous,xi_old,xi_new)
%SMA.UPDATE_POPULATION_ACTIVITY Candidate-coupled threshold traversal.

previous=previous(:);
n=p.Y_activity_population_bins;
if numel(previous)~=n || any(~isfinite(previous)) || any(previous<0)
    error('sma:update_population_activity:InvalidState', ...
        'Previous population activity is invalid.');
end
if p.Y_activity_population_mix<=0
    activity=previous;
    return;
end
q=((1:n).'-0.5)/n;
width=p.Y_activity_population_width;
occupancy_old=0.5*(1.0+tanh((xi_old-q)/(2.0*width)));
occupancy_new=0.5*(1.0+tanh((xi_new-q)/(2.0*width)));
activity=previous+abs(occupancy_new-occupancy_old);
end
