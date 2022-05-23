function Phk_Xk = linear_warp_est(Phprev_Xprev,xk,xprev,Kprev,Ak,Qk,sigk,hrange,w,T,var_scale,area_prec)

% Solves recursive equations to estimate pdf of single parameter (slope) of
% warping function

% INPUTS:
%   Phprev_Xprev - function handle for P(hk-1|Xk-1), or for P(h0) if k=1
%   xk - current state estimate
%   xprev - previous state estimate
%   Kprev - previous state knots, together with xprev define the previous
%       fitted values
%   Ak - state transition matrix
%   Qk - current state noise covariance matrix
%   sigk - variance of P(hk+1|hk)
%   hrange = [hmin hmax] - bounds of slope parameter (default = [0 10])
%   w - function handle for template function
%   T - period of w
%   var_scale - specifies how much to scale variance given by unwarped
%       template residual (i.e. how important is this residual)
%   area_prec - integer value, larger means more accurate rescaling, but
%       more expensive

if isempty(hrange)
    hrange = [0,10];
end
hmin = min(hrange);
hmax = max(hrange);
if isempty(area_prec)
    area_prec = 10;
end

% Define warp transfer pdf, hprev -> hk
Phk_hprev = @(hk,hprev) normpdf(hk,hprev,sigk);

% Define and rescale time update distribution
Phk_Xprev_integrand = @(hk,hprev) (Phk_hprev(hk,hprev)*Phprev_Xprev(hprev));
Phk_Xprev_temp = @(hk) integral(@(hprev) Phk_Xprev_integrand(hk,hprev),hmin,hmax);
Phk_Xprev_area = integral(Phk_Xprev_temp,hmin,hmax,'ArrayValued',true);
Phk_Xprev = @(hk) Phk_Xprev_temp(hk)/Phk_Xprev_area; % nonzero for hk in [hmin,hmax]

% Define and rescale measurement update distribution
% - Certainty (variance) of state given warping and previous state
%       depends on how well previous state matches with unwarped template
% - Assume that a "match" is measured in least squares sense, where there is
%       assumed to be no amplitude variability
Tx = @(h) T/h;
hfunc = @(h,s) h*s;
f_integrand = @(x,K,h,s) (fnval(spmak(K,x'),s) - w(hfunc(h,s)))^2;
varfunc_temp = @(h) integral(@(s) f_integrand(xprev,Kprev,h,s),0,Tx(h),'ArrayValued',true);
varfunc = @(h) var_scale * varfunc_temp(h) / Tx(h); % average residual, rescaled by constant

% Define multivariate probability distributions p(xk|hk,xprev)
Pxk_hkxprev = @(xk,xprev,hk,A,Q) mvnpdf(xk',(A*xprev)',varfunc(hk)*Q);

% Now can define posterior p(hk|Xk)
Phk_Xk_temp = @(h) (Pxk_hkxprev(xk,xprev,h,Ak,Qk)*Phk_Xprev(h)); % Proportional to this

% Find area between [hmin,hmax] to rescale
hvals = linspace(hmin,hmax,area_prec);
postvals = arrayfun(Phk_Xk_temp,hvals);
area = trapz(hvals,postvals);

%ph1x1_area = integral(ph1x1_temp,hmin,hmax,'ArrayValued',true,'RelTol', 1e-4, 'AbsTol', 1e-6); % Assumed to be nonzero for [hmin,hmax]
Phk_Xk = @(h) Phk_Xk_temp(h) / area;


end

