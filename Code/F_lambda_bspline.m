function cost = F_lambda_bspline(knots,coefs,lambda,sig,template,t,T,varargin)

[h,~,w_fun] = h_bspline(knots,coefs,t,T);
t0_template = linspace(0,T,length(template));

w_val = fnval(w_fun,t0_template);
sig_warped = eval_warp(sig,t,h,length(t0_template));


if size(template) ~= size(sig_warped)
    template = template';
end

wts = ones(size(sig_warped));
if nargin > 7
    wts = varargin{1};
    if length(wts) ~= length(sig_warped)
        wts = resample(wts,length(sig_warped),length(wts));
    end
    if size(wts) ~= size(sig_warped)
        wts = wts';
    end
end


% Sum of squared errors term
SSE = sum(wts.*(sig_warped-template).^2);

% Curvature control term
wsq_int = trapz(w_val.^2) * mean(diff(t0_template));

% optional: minimize using unidimensionality

% Compute cost
cost = SSE + lambda * wsq_int;


end