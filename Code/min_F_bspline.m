function [argmin, mincost] = min_F_bspline(knots,coefs_init,lambda,sig,template,t,T,varargin)

if length(coefs_init) >= length(knots)
    error('Must be more knots than coefficients')
end

if nargin > 7
    wts = varargin{1};
    F_lambda = @(coefs) F_lambda_bspline(knots,coefs,lambda,sig,template,t,T,wts);
    argmin = fminsearch(F_lambda,coefs_init);
    mincost = F_lambda(argmin);
else
    F_lambda = @(coefs) F_lambda_bspline(knots,coefs,lambda,sig,template,t,T);
    argmin = fminsearch(F_lambda,coefs_init);
    mincost = F_lambda(argmin);
end
end

