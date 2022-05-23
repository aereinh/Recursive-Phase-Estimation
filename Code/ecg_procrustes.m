function [ecg_frames_reg,t_reg,final_template,H,t_H,costs,prev_template] = ecg_procrustes(ecg_frames,knots,ncoefs,lambda,T,iter,varargin)


% INPUTS
%       ecg_frames - Matrix whose columns contain each frame/beat of
%           ecg_signal, padded with NaN elements to allow for different
%           beat lengths
%       knots - Knot vector for spline basis
%       ncoefs - Number of spline coefficients to be estimated
%       lambda - Regularization parameter
%       T - Period of template function
%       iter - Number of iterations of Procrustes to be performed
%       [optional] wts - Weights for performing weighted least squares 
% 
% OUTPUTS
%       ecg_frames_reg - Registered beat matrix
%       t_reg - Time vector for registered beats
%       final_template - Template calculated with registered beats at
%           current iteration
%       H - Warping function matrix (not normalized)
%       t_H - Time matrix for warping functions
%       costs - Least squares cost between registered curves and template
%       prev_template - Template calculated with registered beats at
%           previous iteration




if ~isempty(knots) && ncoefs >= length(knots)
    error('Number of coefficients must be less than number of knots');
end

% optional parameters
wts = 1;
if nargin > 6
    wts = varargin{1};
    if isempty(wts)
        wts = 1;
    end
end


coefs_init = zeros(1,ncoefs);

ecg_frames_reg = zeros(size(ecg_frames));
t_reg = linspace(0,T,size(ecg_frames,1));
H = nan(size(ecg_frames));
t_H = nan(size(ecg_frames));
costs = nan(1,size(ecg_frames,2));

ecg_frames_rs = nan(size(ecg_frames));
max_len = size(ecg_frames,1);
for k = 1:size(ecg_frames,2)
    frame_k = ecg_frames(:,k);
    frame_k_trunc = frame_k(~isnan(frame_k));
    frame_len = length(frame_k_trunc);
    
    % Resample using splines
    t_or = 1:frame_len;
    t_new = linspace(1,frame_len,length(ecg_frames_rs(:,k)));
    ecg_frames_rs(:,k) = spline(t_or,frame_k_trunc,t_new);
    
%     ecg_frames_rs(:,k) = resample(frame_k_trunc,max_len,frame_len);
end

% First iteration template function is simply the cross-sectional average
template = mean(ecg_frames_rs,2);

if iter == 0
    ecg_frames_reg = ecg_frames_rs;
    final_template = template;
    for k = 1:size(ecg_frames,2)
        frame_k = ecg_frames(:,k);
        frame_k_trunc = frame_k(~isnan(frame_k));
        t_H(1:length(frame_k_trunc),k) = linspace(0,T,length(frame_k_trunc));
        H(1:length(frame_k_trunc),k) = 1:length(frame_k_trunc);
        cost = sum(wts'.*(ecg_frames_reg(:,k)-template).^2);
        costs(1,k) = cost;
    end
    

else
    for k = 1:size(ecg_frames,2)
        frame_k = ecg_frames(:,k);
        frame_k_trunc = frame_k(~isnan(frame_k));
        t = (0:length(frame_k_trunc)-1);
        
        [coefs_min,min_cost] = min_F_bspline(knots,coefs_init,lambda,frame_k_trunc,template,t,T,wts);
        [hmin,t0] = h_bspline(knots,coefs_min,t,T);
                
        H(1:length(hmin),k) = hmin;
        t_H(1:length(t0),k) = t0;
        frame_k_reg = eval_warp(frame_k_trunc,t,hmin,length(template));
        ecg_frames_reg(:,k) = frame_k_reg;
        costs(1,k) = min_cost;
        
    end
    
    prev_template = template;
    new_template = mean(ecg_frames_reg,2);
    completed_iter = 1;

    while completed_iter < iter
        for k = 1:size(ecg_frames,2)
            frame_k = ecg_frames(:,k);
            frame_k_trunc = frame_k(~isnan(frame_k));
            t = (0:length(frame_k_trunc)-1);
            
            [coefs_min,min_cost] = min_F_bspline(knots,coefs_init,lambda,frame_k_trunc,new_template,t,T,wts);
            [hmin,t0] = h_bspline(knots,coefs_min,t,T);
            
            H(1:length(hmin),k) = hmin;
            t_H(1:length(t0),k) = t0;
            frame_k_reg = eval_warp(frame_k_trunc,t,hmin,length(new_template));
            ecg_frames_reg(:,k) = frame_k_reg;
            costs(1,k) = min_cost;
            
        end
        
        prev_template = new_template;
        new_template = mean(ecg_frames_reg,2);
        completed_iter = completed_iter + 1;
    end
    
    final_template = new_template;
end



end


