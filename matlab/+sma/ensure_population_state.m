function state = ensure_population_state(p,state)
%SMA.ENSURE_POPULATION_STATE Add/validate the versioned population state.

n=p.Y_activity_population_bins;
if ~isfield(state,'Y_population_activity')
    if p.Y_activity_population_mix>0
        error('sma:ensure_population_state:MissingActiveState', ...
            ['An enabled threshold-population model requires a v2 state; ' ...
             'migrate the legacy state explicitly before continuing.']);
    end
    state.Y_population_activity=zeros(n,1);
end
activity=state.Y_population_activity(:);
if numel(activity)~=n || any(~isfinite(activity)) || any(activity<0)
    error('sma:ensure_population_state:InvalidState', ...
        'Population activity must contain %d finite nonnegative entries.',n);
end
state.Y_population_activity=activity;
end
