%% Curve registration using ECG Signal

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

figure(1)
clf
plot(ecg_frames)





%% Dynamic Time Warping -- Pairwise

% Take two frames as example
frame1 = ecg_frames(:,1);
frame1 = frame1(~isnan(frame1));

frame2 = ecg_frames(:,2);
frame2 = frame2(~isnan(frame2));

[dist,ix,iy] = dtw(frame1,frame2,10);

% Find midpoint
[xmids,ymids] = midDTW(frame1,frame2,10);

figure(2)
clf
subplot(211)
title('Original signals')
hold on
plot(frame1)
plot(frame2)
subplot(212)
title('DTW')
hold on
plot(frame1(ix))
plot(frame2(iy))

figure(3)
clf
hold on
plot(ix)
plot(iy)

figure(4)
clf
hold on
plot(frame1,'r-')
plot(frame2,'b-')
for i = 1:length(ix)
    plot([ix(i),iy(i)],[frame1(ix(i)),frame2(iy(i))])
end
plot(xmids,ymids,'g-')

%% Dynamic Time Warping -- Multiple signals
% Use this to try and generate template signal

% Start with first two signals/frames
frame1 = ecg_frames(:,1);
frame1 = frame1(~isnan(frame1));

frame2 = ecg_frames(:,3);
frame2 = frame2(~isnan(frame2));

[dist,ix,iy] = dtw(frame1,frame2,10);

mean_pair1 = mean([frame1(ix),frame2(iy)],2);


% Find all pairwise DTW midpoints (templates)
win_len = 10;
npairs = nchoosek(size(ecg_frames,2),2);
pairs = nchoosek(1:size(ecg_frames,2),2);
pw_templates_x = nan(1,npairs);
pw_templates_y = nan(1,npairs);
wts = nan(1,npairs);

for k = 1:npairs
    i = pairs(k,1);
    j = pairs(k,2);
    frame_i = ecg_frames(:,i);
    frame_i = frame_i(~isnan(frame_i));
    frame_j = ecg_frames(:,j);
    frame_j = frame_j(~isnan(frame_j));
    dist_ij = dtw(frame_i,frame_j,win_len);
    [xmids_ij,ymids_ij] = midDTW(frame_i,frame_j,win_len);
    % mean_ij = mean([frame_i(in_i),frame_j(in_j)],2);
    wts(k) = 1/dist_ij;
        
    if length(xmids_ij) > size(pw_templates_x,1)
        new_rows = nan(length(xmids_ij)-size(pw_templates_x,1),npairs);
        pw_templates_x = [pw_templates_x; new_rows];
        pw_templates_y = [pw_templates_y; new_rows];
    end
    % pw_templates(1:length(mean_ij),k) = mean_ij;
    pw_templates_x(1:length(xmids_ij),k) = xmids_ij;
    pw_templates_y(1:length(ymids_ij),k) = ymids_ij;
end

% Resample each midpoint curve to be (1) uniformly spaced and (2) all with
% the same boundary points

pw_templates = nan(size(pw_templates_x));

for k = 1:npairs
    xtemplate_k = pw_templates_x(:,k);
    ytemplate_k = pw_templates_y(:,k);
    template_j = resample(ytemplate_k,xtemplate_k,1,'spline');
    if length(template_j) > size(pw_templates,1)
        new_rows = nan(length(template_j)-size(pw_templates,1),npairs);
        pw_templates = [pw_templates; new_rows];
    end
    pw_templates(1:length(template_j),k) = template_j;
end

pw_templates_rs = nan(size(pw_templates));
max_len = size(pw_templates,1);

for j = 1:npairs
    template_j = pw_templates(:,j);
    template_j_trunc = template_j(~isnan(template_j));
    template_len = length(template_j_trunc);
    pw_templates_rs(:,j) = resample(template_j_trunc,max_len,template_len);
end
    

        
% % Now resample all pairwise templates so that they are same length
% pw_templates_rs = nan(size(pw_templates));
% max_len = size(pw_templates,1);
% for k = 1:npairs
%     template_k = pw_templates(:,k);
%     template_k_trunc = template_k(~isnan(template_k));
%     template_len = length(template_k_trunc);
%     pw_templates_rs(:,k) = resample(template_k_trunc,max_len,template_len);
% end
% 

% Find mean of all midpoint DTW templates
overall_template = mean(pw_templates_rs,2);

% Compute a weighted mean using DTW distances (adjust so weights add up to
% n=120)
wts_adj = (npairs/sum(wts)) * wts;

wt_templates = (pw_templates_rs'.*wts_adj')';

overall_template_wt = mean(wt_templates,2);

figure(1)
clf
hold on
p1 = plot(pw_templates_rs,'LineWidth',.5);

p2 = plot(overall_template,'r-','LineWidth',2);
p2.Color(4) = 1;
p3 = plot(overall_template_wt,'y-','LineWidth',2);
p3.Color(4) = 1;
    


%% Cruder version: No DTW, but just resampling
ecg_frames_rs = nan(size(ecg_frames));
max_len = size(ecg_frames,1);
for k = 1:size(ecg_frames,2)
    frame_k = ecg_frames(:,k);
    frame_k_trunc = frame_k(~isnan(frame_k));
    frame_len = length(frame_k_trunc);
    ecg_frames_rs(:,k) = resample(frame_k_trunc,max_len,frame_len);
end

overall_template2 = mean(ecg_frames_rs,2);

figure(2)
clf
hold on
plot(ecg_frames_rs,'LineWidth',1)
plot(overall_template2,'r-','LineWidth',2)


%% FDA using Procrustes algorithm and B-Spline bases

% Start off by resampling all frames to get same length
ecg_frames_rs = nan(size(ecg_frames));
max_len = size(ecg_frames,1);
for k = 1:size(ecg_frames,2)
    frame_k = ecg_frames(:,k);
    frame_k_trunc = frame_k(~isnan(frame_k));
    frame_len = length(frame_k_trunc);
    ecg_frames_rs(:,k) = resample(frame_k_trunc,max_len,frame_len);
end

% First iteration template function is simply the cross-sectional average
template1 = mean(ecg_frames_rs,2);

figure(1)
clf
hold on
plot(ecg_frames_rs,'--','LineWidth',.5)
plot(template1,'r','LineWidth',2)


% Now for each frame i, we want to find optimal wi(t) with warping function
% hi(t)=c0+c1*integral{exp(wi(u))du} so that
% xi{hi(t)} closely resembles initial template y(t)

fs = 216;

% Take frame1 for example
frame1 = ecg_frames(~isnan(ecg_frames(:,1)),1)';
t1 = (0:length(frame1)-1);

% Specify knots for B-spline
knots = [1, 13, 20, 100, 162, 260, 309, 337, 347, 359]/359;
coefs_init = [0,0,0,0,0,0,0];

[h_init,t0,w] = h_bspline(knots,coefs_init,t1,1);
t0_template = linspace(0,1,length(template1));
frame1_nowarp = eval_warp(frame1,h_init,359);
cost_init = F_lambda_bspline(knots,coefs_init,.01,frame1,template1,t1,1);

% Now register curve by minimizing squared error and curvature terms
lambda = 0.01;
[coefs_min,cost_min] = min_F_bspline(knots,coefs_init,lambda,frame1,template1,t1,1);
hmin = h_bspline(knots,coefs_min,t1,1);
frame1_reg = eval_warp(frame1,hmin,359);

figure(3)
clf
subplot(211)
hold on
% plot(t0,frame1,'b')
% plot(t0_template,frame1_randwarp,'r-')
plot(t0_template,frame1_reg,'r*')
plot(t0_template,frame1_nowarp,'b--','LineWidth',2)
plot(t0_template,template1,'black','LineWidth',2)

subplot(212)
hold on
plot(t0,hmin,'r*')
plot(t0,h_init,'b--','LineWidth',2)


%% Use Procrustes algorithm to obtain registered curves and templates
test_frames = ecg_frames;

[~,t_reg,template0] = ecg_procrustes(test_frames,[],[],[],1,0);

% Uniform knots
nknots = 5;
temp = linspace(t_reg(1),t_reg(end),nknots+2);
knots = temp(2:end-1);

% Uniform weights
wts = ones(1,size(ecg_frames,1));
% wts(1:10) = 0.5; wts(end-10:end) = 0.5;

ncoefs = nknots - 3;
lambda = 1;
T = 1;

% Resampling only
[ecg_frames_reg0,t_reg0,final_template0,H0,t_H0,costs0] = ecg_procrustes(test_frames,knots,ncoefs,lambda,T,0,wts);

% One iteration of Procrustes
[ecg_frames_reg1,t_reg1,final_template1,H1,t_H1,costs1] = ecg_procrustes(test_frames,knots,ncoefs,lambda,T,1,wts);

% Two iterations of Procrustes
%[ecg_frames_reg2,t_reg2,final_template2,H2,t_H2,costs2] = ecg_procrustes(test_frames,knots,ncoefs,lambda,T,2,wts);


figure(1)
clf
subplot(311)
hold on
title('iter=0, lambda = 0.0001')
plot(t_reg0,ecg_frames_reg0,'--','LineWidth',.5)
plot(t_reg0, final_template0,'r','LineWidth',1)
subplot(312)
hold on
title('iter=1')
plot(t_reg1,ecg_frames_reg1,'--','LineWidth',.5)
plot(t_reg1, final_template1,'r','LineWidth',1)
subplot(313)
hold on
title('iter=2')
plot(t_reg2,ecg_frames_reg2,'--','LineWidth',.5)
plot(t_reg2, final_template2,'r','LineWidth',1)

figure(2)
clf
subplot(311)
hold on
title('iter=0,lambda=0.0001')
plot(t_H0,H0)
subplot(312)
hold on
title('iter=1')
plot(t_H1,H1)
subplot(313)
hold on
title('iter=2')
plot(t_H2,H2)

figure(3)
clf
subplot(211)
hold on
plot(t_reg2,final_template2,'black')
plot(t_reg0,ecg_frames_reg0(:,3),'r--')
%plot(t_reg1,ecg_frames_reg1(:,3),'g--')
plot(t_reg2,ecg_frames_reg2(:,3),'b--')
legend('template','resampled curve','registered curve')
subplot(212)
hold on
plot(t_H0(:,3),H0(:,3),'r')
plot(t_H2(:,3),H2(:,3),'b')
legend('resampled (linear)','registered (nonlinear)')



%% Curve Registration with Endpoint Constrained EKF


figure(1)
clf
plot(t_reg0, final_template0)

f = @(t) spline(t_reg0,final_template0,mod(t,1));
df = @(t) fnval(fnder(spline(t_reg0,final_template0)), mod(t,1));
t = linspace(0,5,1000);

% figure(2)
% clf
% hold on
% plot(t,f(t))
% plot(t,df(t))

d = 3;
J = 4;
x0 = zeros(J+2,1); x0(J+2) = 1;
P0 = diag([repmat(1e2,1,J) 1e1 1e-12]);

i = 2;

xdat = test_frames(~isnan(test_frames(:,i)),i);
tdat = linspace(0,1,length(xdat));
ndat = length(tdat);
Q = diag(1e-12*ones(1,J+2));
R = 1;
ytild = f;
dytild = df;
[hinv,~,x1_n] = phase_function_estimator_constr(x0,P0,tdat,xdat,d,J,Q,R,ytild,dytild,ndat,2,0,inf,0);
xfit = x1_n(J+1,end) + x1_n(J+2,end) * ytild(hinv(tdat));


% curve fit residuals



figure(1)
clf
subplot(211)
hold on
plot(tdat,xdat,'k.')
plot(tdat,xfit)
plot(tdat,ytild(tdat),'--')
subplot(212)
plot(tdat,hinv(tdat))



% figure(1)
% clf
% plot(test_frames)



%% Functions


function [phase_func,freq_func,x1_n,P1_n,K,zfit,phase1_n] = phase_function_estimator_constr(x0,P0,tdat,xdat,d,J,Q,R,ytild,dytild,varargin)

Ktemp = -d:J;
K = (tdat(end)+.0001)*Ktemp/Ktemp(end-d);
K(1:d+1) = linspace(K(d+1)-.01,K(d+1),d+1);
K(J+1:end) = linspace(K(J+1),K(J+1)+.01,d+1);

xprev = x0;
Pprev = P0;

ndat = length(tdat);
if nargin > 10
    ndat = varargin{1};
end

zfit = zeros(1,ndat);
phase1_n = zeros(1,ndat);
x1_n = zeros(J+2,ndat);
P1_n = [];

for k = 1:ndat
    [xk,Pk,~,~,~,~,~,~,~,Ck] = kalman_warp4(K,xprev,Pprev,R,tdat(k),xdat(k),Q,ytild,dytild);
    
    % increasing coefficients constraint
    if nargin > 11 & varargin{2} == 1
        for i = 2:J
            if xk(i) < 0
                xk(i) = 0;
            end
        end 
    end

    % derivative constraint
    if nargin > 11 & varargin{2} == 2
        lb = 0;
        ub = inf;
        if nargin > 12
            lb = varargin{3};
        end
        if nargin > 13
            ub = varargin{4};
        end
            
        
        muk = find(K <= tdat(k),1,'last');
        bspline = d * make_B_prime(muk,d,K);
        for delta = d-1:-1:1
            bspline = make_B(muk,delta,tdat(k),K) * bspline;
        end
        Cderiv = [zeros(1,muk-(d+1)) bspline zeros(1,(J-muk))];
        D_c = flip(cumsum(flip(Cderiv)));
        %W = 1;
        W = inv(Pk(1:J,1:J));
        xk_c = xk(1:J);
        if D_c*xk_c <= lb
            xk_c = xk(1:J) - inv(W) * D_c' * inv(D_c*inv(W)*D_c')*(D_c*xk(1:J)-lb);
        elseif D_c*xk_c >= ub
            xk_c = xk(1:J) - inv(W) * D_c' * inv(D_c*inv(W)*D_c')*(D_c*xk(1:J)-ub);
        end
        xk(1:J) = xk_c;
    end
    
    % endpoint constraints
    if nargin > 14
        C = make_C(K,tdat(1),d,J);
        D_c = flip(cumsum(flip(C)));
        d_c = varargin{5};
        if nargin > 15
            C2 = make_C(K,tdat(end),d,J);
            D_c = [flip(cumsum(flip(C))); flip(cumsum(flip(C2)))];
            d_c = [varargin{5};varargin{6}];
        end
        W = inv(Pk(1:J,1:J));
        xk_c = xk(1:J) - inv(W) * D_c' * inv(D_c*inv(W)*D_c')*(D_c*xk(1:J)-d_c);
        xk(1:J) = xk_c;
    end
    

    x1_n(:,k) = xk;
    P1_n(:,:,k) = Pk;
    
    hk = cumsum(xk(1:J));
    
    zfit(k) = xk(J+1) + xk(J+2) * ytild(Ck*hk);
    phase1_n(k) = Ck*hk;
    
    xprev = xk;
    Pprev = Pk;
end

%hn_func = @(t) fnval(spmak(K,hk(1:J)'),t);
phase_func = @(t) fnval(spmak(K,hk(1:J)'),t);
freq_func = @(t) fnval( fnder( spmak(K,hk(1:J)')),t);

end
