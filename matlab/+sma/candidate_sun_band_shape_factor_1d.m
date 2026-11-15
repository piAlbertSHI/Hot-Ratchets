function [H,meta]=candidate_sun_band_shape_factor_1d(length_ratio,nu)
%SMA.CANDIDATE_SUN_BAND_SHAPE_FACTOR_1D
% Evaluate Sun/Zhong (2000) Eq.(28)--(30) for the elastic cylindrical
% inclusion shape factor.  This is a pure validation kernel; it contains no
% phase kinetics, thermal law, observation access or state mutation.

assert(isnumeric(length_ratio) && isscalar(length_ratio) && isfinite(length_ratio) && ...
    length_ratio>0,'sma:sun_shape:Domain','length_ratio must be positive.');
assert(isnumeric(nu) && isscalar(nu) && isfinite(nu) && nu>-1 && nu<0.5, ...
    'sma:sun_shape:Domain','nu must be in (-1,0.5).');
lambda=double(length_ratio); n=double(nu);
% The infinite-interval adaptive splitter can exhaust its interval budget on
% the oscillatory tail.  Use deterministic composite trapezoids at two long
% cutoffs and expose their difference as a tail-convergence diagnostic.
step=0.05;
t4000=(0:step:4000).'; t2000=(0:step:2000).';
I2000=trapz(t2000,shape_integrand(t2000,lambda,n));
I4000=trapz(t4000,shape_integrand(t4000,lambda,n));
I=I4000; err=abs(I4000-I2000); H=1+I;
meta=struct('source_equations','SunZhong2000_Eq28_Eq30', ...
    'length_ratio',lambda,'nu',n,'integral_value',I, ...
    'quadrature_error_estimate',err, ...
    'tail_convergence_error_2000_to_4000',err, ...
    'integration_step',step,'integration_cutoffs',[2000 4000], ...
    'validation_only',true, ...
    'source_authoritative',false,'production_connected',false, ...
    'state_modified',false,'observation_accessed',false);
assert(isfinite(H) && H>0,'sma:sun_shape:Numerics', ...
    'Sun/Zhong shape factor is nonpositive or nonfinite.');
end

function y=shape_integrand(t,lambda,nu)
t=double(t); y=zeros(size(t));
small=(t<1e-6);
% The removable t=0 limit of G(t)sin^2(lambda*t/2) is -lambda/pi.
y(small)=-lambda/pi;
idx=~small;
if any(idx(:))
    x=t(idx);
    % Scaled modified Bessel functions avoid overflow in I0/I1 at infinity.
    ratio=besseli(0,x,1)./besseli(1,x,1);
    rho=2*(1-nu)+x.*ratio;
    denominator=(1-2*nu-rho).*ratio+x+rho./x;
    G=8*(1+nu)./(pi*lambda*x.^3)./denominator;
    y(idx)=G.*sin(x*lambda/2).^2;
end
end
