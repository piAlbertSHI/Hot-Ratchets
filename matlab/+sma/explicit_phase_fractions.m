function [xi_A,xi_R,xi_M,relative_R] = explicit_phase_fractions( ...
    p,xi_M,seq,T,training_activity)
%SMA.EXPLICIT_PHASE_FRACTIONS Zhang Eqs. (7): A/R/M allocation.
%
% The existing R law supplies the R fraction inside the currently
% transformable A+R inventory.  The explicit candidate converts that
% relative fraction to an absolute fraction and enforces
% xi_A+xi_R+xi_M=1 exactly.

if nargin<5
    training_activity=0.0;
end
xi_M=min(max(xi_M,0.0),1.0);
relative_R=sma.r_phase_fraction(p,seq,T,training_activity);
relative_R=min(max(relative_R,0.0),1.0);
xi_R=(1.0-xi_M)*relative_R;
xi_A=1.0-xi_M-xi_R;
xi_A=max(xi_A,0.0);
assert(abs(xi_A+xi_R+xi_M-1.0)<5.0e-14 && ...
    min([xi_A,xi_R,xi_M])>=-1.0e-14, ...
    'sma:explicit_phase_fractions:Conservation', ...
    'Explicit A/R/M phase allocation is not conservative.');
end
