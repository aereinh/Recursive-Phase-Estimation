%% Read in data

load('ECGSignal')
ecg = ECGSignal(:,1);
L = length(ecg);
t = 0:100/L:100-100/L;
fs = 1/mean(diff(t));

L_samp = 5000;
ecg_samp = ecg(1:L_samp);
t_samp = t(1:L_samp);
ecg_filt = lowpass(ecg_samp,10,fs);

[pks,locs] = findpeaks(ecg_filt,'MinPeakHeight',.2);
max_interval = max(diff(locs));

ecg_frames = NaN(max_interval,length(pks)-1);

for i = 1:length(pks)-1
    frame_len = locs(i+1)-locs(i)+1;
    ecg_frames(1:frame_len,i) = ecg_filt(locs(i):locs(i+1));
end

ecg_frames(end,:) = NaN;

%% Fitting B-Splines to data

test_ecg = ecg_frames(~isnan(ecg_frames(:,1)),1);
L = length(test_ecg);
test_t = t(1:L);

% Fit a B-Spline using global criteria
deg = 3;
nknots = 5;
temp_knots = linspace(test_t(1),test_t(end),nknots+2);
knots = temp_knots(2:end-1);
ncoefs = nknots-deg;

coefs_init = zeros(1,ncoefs);

%% 

[spl_vals,min_cost,opt_coefs] = fit_bspline(test_ecg',test_t,knots,coefs_init);

function cost = bspline_cost(x,t,knots,coefs)
    bspl = spmak(knots,coefs);
    spl_val = fnval(bspl,t);
    cost = sum((spl_val - x).^2);
end

function [spl_vals, min_cost, opt_coefs] = fit_bspline(x,t,knots,coefs_init)
    cost_fun = @(coefs) bspline_cost(x,t,knots,@coefs);
    [min_cost,opt_coefs] = fminsearch(cost_fun,coefs_init);
    spl = spmak(knots,opt_coefs);
    spl_vals = fnval(spl,t);
end


