function [h, t0, w] = h_bspline(knots,coefs,t,T0)

% INPUTS:
%   knots -- sequence of breakpoints for B-spline
%   coefs -- coefficients for B-spline
%   t -- time vector for which B-spline is defined
%   T0 -- specified end time such that h(t(end)) = T0

% Define w(t) to be B-spline with given knots and coefficients
w = spmak(knots,coefs);
Ti = t(end)-t(1);

% Find integral of that B-spline (call this W(t))
w_int = fnint(w);

% Find numerical integral of exp(W(t)) over 0 to T0, use this to calculate
% constant c1
nt = 10000;
t_temp = linspace(0,T0,nt);
w_int_vals = exp(spval(w_int,t_temp));
c1 = Ti / (sum(w_int_vals)*mean(diff(t_temp)));

% With c1, we can find values of h(t) for t0 in [0,T0]
t0 = linspace(0,T0,length(t));
h = c1 * cumsum(exp(spval(w_int,t0))) * mean(diff(t0));
h = h - h(1);

end