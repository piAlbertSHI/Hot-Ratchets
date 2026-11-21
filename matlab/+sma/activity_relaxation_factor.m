function factor = activity_relaxation_factor(p,activity,peak_fraction)
%SMA.ACTIVITY_RELAXATION_FACTOR Continuous fast/slow defect relaxation.
% Every material uses the same equation.  The legacy single exponential is
% recovered exactly when fast_weight=1 and fast_rate=lambda_Y.

activity=max(activity,0.0);
if nargin<3, peak_fraction=0.0; end
peak_fraction=min(max(peak_fraction,0.0),1.0);
w=p.Y_activity_fast_weight;
legacy=exp(-p.lambda_Y*activity);
spectrum=w*exp(-p.Y_activity_fast_rate*activity)+ ...
    (1.0-w)*exp(-p.Y_activity_slow_rate*activity);
% Partial transformation samples only part of the resistance population.
% This continuous capacity weighting contains no material or protocol switch.
partial_weight=(1.0-peak_fraction)^p.Y_activity_partial_power;
factor=(1.0-partial_weight)*legacy+partial_weight*spectrum;
end
