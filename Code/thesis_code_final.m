% Simulation 1: Linear phase Sine wave
T = 2*pi; % period
ytild = @(t) sin(T*t).^1; % template curve
dytild = @(t) 1*T*sin(T*t).^0*cos(T*t);

Tx = 1;

htrue = @(t) 10*t-.2;%2*(t+2*t.^2);%5*t.^2 + t + .1;%5*t.^2; % true inverse warping/phase function
atrue = [0 1]; % true amplitude parameters
anoise = .1; % noise amplitude
L = 5000;

rng('default')
rng(999)
tdat = linspace(0,Tx,L)';
xdat = atrue(1)+atrue(2)*ytild(htrue(tdat))+anoise*randn(L,1);
fs = 1./mean(diff(tdat));
dt = mean(diff(tdat));

% Method 1: Hilbert transform estimate
hilb = hilbert(xdat);
ph_hilb_wr = angle(hilb);
ph_hilb_unwr = unwrap2(ph_hilb_wr',T,.5*T)';
ph_hilb_norm = (ph_hilb_unwr + pi/2)/(2*pi);
hilb_ph_recon = cos(ph_hilb_unwr);

% Method 2: DFT (works for linear)
% xfft = fft(xdat);
% fz = linspace(0,fs,length(xdat));
% Lwin = round(length(xdat)/10);
% beta = 3;
% win = kaiser(Lwin,beta);
% noverlap = 0;%Lwin - 1;
% frange = [0 30];
% [z,fz,tz] = czt_spectrogram(xdat-mean(xdat),kaiser(Lwin,beta),noverlap,[],frange,fs);
% z_adj = adj_phase(z,tz,fz);
% pkfreqs = zeros(1,size(z,2));
% phoffsets1 = zeros(1,size(z,2));
% phoffsets2 = zeros(1,size(z,2));
% for i = 1:size(z,2)
%     [pk,pkfreq_loc] = max(abs(z(:,i)));
%     pkfreqs(i) = fz(pkfreq_loc);
%     phoffsets1(i) = angle(z(pkfreq_loc,i));
%     phoffsets2(i) = angle(z_adj(pkfreq_loc,i));
% end

% EKF 1: Raw Phase Estimator
x0 = zeros(4,1); x0(3) = 1;
P0 = diag([1e2 1e2 1e-12 1e-12]);
Q = diag([1e-12 1e-12 1e-12 1e-12]);
R = 1;
x1_n_rawekf = raw_phase_estimator(x0,P0,xdat,tdat,Q,R,ytild,dytild,0);


% EKF 2: Phase Function Estimator
d = 1;
J = 10;
x0 = zeros(J+2,1);
x0(J+2) = 1;
P0 = diag([1e2 repmat(1e2,1,J-1) 1e1 1e1]);
Q = diag([1e-12 repmat(1e-12,1,J-1) 1e-12 1e-12]); 
R = 1;
ndat = length(tdat);
[hn_func,~,x1_n_func,~,K,zfit,phase1_n_ekf1] = phase_function_estimator_constr(x0,P0,tdat,xdat,d,J,Q,R,ytild,dytild,round(ndat*1),2,2,15);
yn_func = @(t) x1_n_func(end,end)*ytild(hn_func(t))+x1_n_func(end-1,end);

figure(1)
clf
subplot(311)
hold on
plot(tdat,htrue(tdat),'black')
plot(tdat,ph_hilb_norm,'r')
plot(tdat,x1_n_rawekf(1,:),'b')
plot(tdat,hn_func(tdat),'b')
%plot(tdat,zfit,'g')
%legend('true','hilbert','ekf1','ekf2')
title('Phase Estimates')
subplot(312)
hold on
plot(tdat,ph_hilb_norm-htrue(tdat),'r')
plot(tdat,x1_n_rawekf(1,:)'-htrue(tdat),'b')
plot(tdat,hn_func(tdat)-htrue(tdat),'g')
yline(0,'--')
%legend('hilbert','ekf1','ekf2')
title('Phase Errors')
subplot(313)
hold on
plot(tdat,xdat,'r.','MarkerSize',.1)
%plot(tdat,hilb_ph_recon,'r')
% plot(tdat,fourier_ph_recon,'b')
%plot(tdat,ytild(x1_n_rawekf(1,:)),'b','LineWidth',1)
plot(tdat,yn_func(tdat),'g','LineWidth',2)
xline(tdat(ndat))
%legend('data','hilbert','ekf1','ekf2')
title('Phase Reconstructed Signals')

mse_hilb = sum( (ph_hilb_norm-htrue(tdat)).^2) / length(tdat);
mse_ekf1 = sum ( (x1_n_rawekf(1,:)'-htrue(tdat)).^2) / length(tdat);
mse_ekf2 = sum ( ( hn_func(tdat) - htrue(tdat) ).^2) / length(tdat);

%% Simulation 2: Quadratic Sine wave
T = 2*pi; % period
ytild = @(t) sin(T*t).^1; % template curve
dytild = @(t) 1*T*sin(T*t).^0*cos(T*t);

Tx = 1;

htrue = @(t) 5*t.^2+t+.1;%2*(t+2*t.^2);%5*t.^2 + t + .1;%5*t.^2; % true inverse warping/phase function
atrue = [0 1]; % true amplitude parameters
anoise = .1; % noise amplitude
L = 5000;

rng('default')
rng(999)
tdat = linspace(0,Tx,L)';
xdat = atrue(1)+atrue(2)*ytild(htrue(tdat))+anoise*randn(L,1);
fs = 1./mean(diff(tdat));
dt = mean(diff(tdat));

% Method 1: Hilbert transform estimate
hilb = hilbert(xdat);
ph_hilb_wr = angle(hilb);
ph_hilb_unwr = unwrap2(ph_hilb_wr',T,.5*T)';
ph_hilb_norm = (ph_hilb_unwr + pi/2)/(2*pi);
hilb_ph_recon = cos(ph_hilb_unwr);

% EKF 1: Raw Phase Estimator
x0 = zeros(4,1); x0(3) = 1;
P0 = diag([1e2 1e2 1e-12 1e-12]);
Q = diag([1e-12 1e0 1e-12 1e-12]);
R = 1;
x1_n_rawekf = raw_phase_estimator(x0,P0,xdat,tdat,Q,R,ytild,dytild,0);


% EKF 2: Phase Function Estimator
d = 1;
J = 10;
x0 = zeros(J+2,1);
x0(J+2) = 1;
P0 = diag([1e2 repmat(1e2,1,J-1) 1e-12 1e-12]);
Q = diag([1e-12 repmat(1e-12,1,J-1) 1e-12 1e-12]); 
R = 1;
ndat = length(tdat);
[hn_func,~,x1_n_func,~,K,zfit,phase1_n_ekf1] = phase_function_estimator_constr(x0,P0,tdat,xdat,d,J,Q,R,ytild,dytild,round(ndat*1),2,2,15);
yn_func = @(t) x1_n_func(end,end)*ytild(hn_func(t))+x1_n_func(end-1,end);

figure(1)
clf
subplot(311)
hold on
plot(tdat,htrue(tdat),'black')
plot(tdat,ph_hilb_norm,'r')
plot(tdat,x1_n_rawekf(1,:),'b')
plot(tdat,hn_func(tdat),'b')
%plot(tdat,zfit,'g')
%legend('true','hilbert','ekf1','ekf2')
title('Phase Estimates')
subplot(312)
hold on
plot(tdat,ph_hilb_norm-htrue(tdat),'r')
plot(tdat,x1_n_rawekf(1,:)'-htrue(tdat),'b')
plot(tdat,hn_func(tdat)-htrue(tdat),'g')
yline(0,'--')
%legend('hilbert','ekf1','ekf2')
title('Phase Errors')
subplot(313)
hold on
plot(tdat,xdat,'r.','MarkerSize',.1)
%plot(tdat,hilb_ph_recon,'r')
% plot(tdat,fourier_ph_recon,'b')
%plot(tdat,ytild(x1_n_rawekf(1,:)),'b','LineWidth',1)
plot(tdat,yn_func(tdat),'g','LineWidth',2)
xline(tdat(ndat))
%legend('data','hilbert','ekf1','ekf2')
title('Phase Reconstructed Signals')

mse_hilb = sum( (ph_hilb_norm-htrue(tdat)).^2) / length(tdat);
mse_ekf1 = sum ( (x1_n_rawekf(1,:)'-htrue(tdat)).^2) / length(tdat);
mse_ekf2 = sum ( ( hn_func(tdat) - htrue(tdat) ).^2) / length(tdat);


%% Simulation 3: B-Spline Phase Sinusoid

T = 2*pi; % period
ytild = @(t) sin(T*t).^1; % template curve
dytild = @(t) 1*T*sin(T*t).^0*cos(T*t);

Tx = 1;

d = 3;
J = 10;
Ktemp = -d:J;
K = (1.0001)* Ktemp / Ktemp(end-d);
x = [-1 0 1 1.2 3.8 4.5 5 7 9 10];
spl = spmak(K,x);
htrue = @(t) fnval(spl,t);

atrue = [0 1]; % true amplitude parameters
anoise = .1; % noise amplitude
L = 5000;

rng('default')
rng(999)
tdat = linspace(0,Tx,L)';
xdat = atrue(1)+atrue(2)*ytild(htrue(tdat))+anoise*randn(L,1);
fs = 1./mean(diff(tdat));
dt = mean(diff(tdat));

% Method 1: Hilbert transform estimate
hilb = hilbert(xdat);
ph_hilb_wr = angle(hilb);
ph_hilb_unwr = unwrap2(ph_hilb_wr',T,.5*T)';
ph_hilb_norm = (ph_hilb_unwr + pi/2)/(2*pi);
hilb_ph_recon = cos(ph_hilb_unwr);

% EKF 1: Raw Phase Estimator
x0 = zeros(4,1); x0(3) = 1;
P0 = diag([1e2 1e2 1e-12 1e-12]);
Q = diag([1e-12 1e0 1e-12 1e-12]);
R = 1;
x1_n_rawekf = raw_phase_estimator(x0,P0,xdat,tdat,Q,R,ytild,dytild,0);

% EKF 2: Phase Function Estimator
d = 3;
J = 10;
x0 = zeros(J+2,1);
x0(J+2) = 1;
P0 = diag([1e2 repmat(1e2,1,J-1) 1e-12 1e-12]);
Q = diag([1e-12 repmat(1e-12,1,J-1) 1e-12 1e-12]); 
R = 1;
ndat = length(tdat);
[hn_func,~,x1_n_func,~,K,zfit,phase1_n_ekf1] = phase_function_estimator_constr(x0,P0,tdat,xdat,d,J,Q,R,ytild,dytild,round(ndat*1),2,2,15);
yn_func = @(t) x1_n_func(end,end)*ytild(hn_func(t))+x1_n_func(end-1,end);

figure(1)
clf
subplot(311)
hold on
plot(tdat,htrue(tdat),'black')
plot(tdat,ph_hilb_norm,'r')
plot(tdat,x1_n_rawekf(1,:),'b')
plot(tdat,hn_func(tdat),'b')
%plot(tdat,zfit,'g')
%legend('true','hilbert','ekf1','ekf2')
title('Phase Estimates')
subplot(312)
hold on
plot(tdat,ph_hilb_norm-htrue(tdat),'r')
plot(tdat,x1_n_rawekf(1,:)'-htrue(tdat),'b')
plot(tdat,hn_func(tdat)-htrue(tdat),'g')
yline(0,'--')
%legend('hilbert','ekf1','ekf2')
title('Phase Errors')
subplot(313)
hold on
plot(tdat,xdat,'r.','MarkerSize',.1)
%plot(tdat,hilb_ph_recon,'r')
% plot(tdat,fourier_ph_recon,'b')
%plot(tdat,ytild(x1_n_rawekf(1,:)),'b','LineWidth',1)
plot(tdat,yn_func(tdat),'g','LineWidth',2)
xline(tdat(ndat))
%legend('data','hilbert','ekf1','ekf2')
title('Phase Reconstructed Signals')

mse_hilb = sum( (ph_hilb_norm-htrue(tdat)).^2) / length(tdat);
mse_ekf1 = sum ( (x1_n_rawekf(1,:)'-htrue(tdat)).^2) / length(tdat);
mse_ekf2 = sum ( ( hn_func(tdat) - htrue(tdat) ).^2) / length(tdat);

%% Simulation 4: Double beat sinusoid

T = 2*pi; % period
ytild = @(t) sin(4*pi*t)+sin(2*pi*t + pi/2);
dytild = @(t) 4*pi*cos(4*pi*t) + 2*pi*cos(2*pi*t+pi/2);

Tx = 1;

d = 3;
J = 10;
Ktemp = -d:J;
K = (1.0001)* Ktemp / Ktemp(end-d);
x = [-1 0 1 1.5 2 2.5 3 3.5 3.5 4.5];
% 0*[0 1 2 2.5 3 3.5 4 4.5 4.5 5.5]
spl = spmak(K,x);
htrue = @(t) fnval(spl,t);

atrue = [0 1]; % true amplitude parameters
anoise = .1; % noise amplitude
L = 5000;

rng('default')
rng(999)
tdat = linspace(0,Tx,L)';
xdat = atrue(1)+atrue(2)*ytild(htrue(tdat))+anoise*randn(L,1);
fs = 1./mean(diff(tdat));
dt = mean(diff(tdat));

% Method 1: Hilbert transform estimate
hilb = hilbert(xdat);
ph_hilb_wr = angle(hilb);
ph_hilb_unwr = unwrap2(ph_hilb_wr',T,.5*T)';
ph_hilb_norm = (ph_hilb_unwr + pi/2)/(2*pi);
hilb_ph_recon = cos(ph_hilb_unwr);

% EKF 1: Raw Phase Estimator
x0 = zeros(4,1); x0(3) = 1;
P0 = diag([1e2 1e2 1e-12 1e-12]);
Q = diag([1e-5 1e0 1e-12 1e-12]);
R = 1;
x1_n_rawekf = raw_phase_estimator(x0,P0,xdat,tdat,Q,R,ytild,dytild,0);

% EKF 2: Phase Function Estimator
d = 3;
J = 10;
x0 = zeros(J+2,1);
x0(J+2) = 1;
P0 = diag([1e2 repmat(1e2,1,J-1) 1e-12 1e-12]);
Q = diag([1e-12 repmat(1e-12,1,J-1) 1e-12 1e-12]); 
R = 1;
ndat = length(tdat);
[hn_func,~,x1_n_func,~,K,zfit,phase1_n_ekf1] = phase_function_estimator_constr(x0,P0,tdat,xdat,d,J,Q,R,ytild,dytild,round(ndat*1),2,1,15);
yn_func = @(t) x1_n_func(end,end)*ytild(hn_func(t))+x1_n_func(end-1,end);

figure(1)
clf
subplot(311)
hold on
plot(tdat,htrue(tdat),'black')
plot(tdat,ph_hilb_norm,'r')
plot(tdat,x1_n_rawekf(1,:),'b')
plot(tdat,hn_func(tdat),'g')
%plot(tdat,zfit,'g')
%legend('true','hilbert','ekf1','ekf2')
title('Phase Estimates')
subplot(312)
hold on
%plot(tdat,ph_hilb_norm-htrue(tdat),'r')
plot(tdat,x1_n_rawekf(1,:)'-htrue(tdat),'b')
plot(tdat,hn_func(tdat)-htrue(tdat),'g')
yline(0,'--')
%legend('hilbert','ekf1','ekf2')
title('Phase Errors')
subplot(313)
hold on
plot(tdat,xdat,'r.','MarkerSize',.1)
%plot(tdat,hilb_ph_recon,'r')
% plot(tdat,fourier_ph_recon,'b')
%plot(tdat,ytild(x1_n_rawekf(1,:)),'b','LineWidth',1)
plot(tdat,yn_func(tdat),'g','LineWidth',2)
xline(tdat(ndat))
%legend('data','hilbert','ekf1','ekf2')
title('Phase Reconstructed Signals')

mse_hilb = sum( (ph_hilb_norm-htrue(tdat)).^2) / length(tdat);
mse_ekf1 = sum ( (x1_n_rawekf(1,:)'-htrue(tdat)).^2) / length(tdat);
mse_ekf2 = sum ( ( hn_func(tdat) - htrue(tdat) ).^2) / length(tdat);



%% Simulation 5 (see thesis_code_snippets2 (recursive curve registration))

%% Simulation 6 (see thesis_code_snippets2 (simulated ecg curve))

%% Simulation 7-8




%% Simulation 9



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
    u = 0; B = 0; A = diag(ones(1,J+2));
    muk = find(K <= tdat(k),1,'last');
    Ck = make_C(K,tdat(k),d,J);
    xpredk = A*xprev;  
    innovk = xdat(k) - (xpredk(J+2) * ytild(Ck*xpredk(1:J)) + xpredk(J+1));
    Hk = [xpredk(J+2)*dytild(Ck*xpredk(1:J))*Ck, 1, ytild(Ck*xpredk(1:J))];
    
    [xk,Pk] = kalman(xprev,Pprev,u,xdat(k),A,B,Hk,Q,R,innovk);
        
    % increasing coefficients constraint
    if nargin > 11 & varargin{2} == 1
        if d > 1
            error('only applies for first-degree B-splines')
        end
        dbar = 0;
        if nargin > 12
            dbar = varargin{3} * max(diff(K));
        end
 
        Ltemp = tril(ones(J));
        L = inv(Ltemp);
        D_c = [L(2:end,:) zeros(J-1,2)];
        d_c = dbar*ones(J-1,1);
        C = diag(ones(1,J+2));
        xk_uc = xk;
        if any(D_c * xk_uc < d_c)
            xk = lsqlin(C,xk_uc,-D_c,-d_c);
        end
        
    end

    % positive amplitude
    if xk(end) < 0
        xk(end) = 0;
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
            
        
        Cderiv = make_C(K,tdat(k),d,J,1);
        D_c = [Cderiv 0 0];
        W = inv(Pk);
        
        xk_c = xk;
        if D_c*xk <= lb
            xk_c = xk - inv(W) * D_c' * inv(D_c*inv(W)*D_c')*(D_c*xk-lb);
        elseif D_c*xk >= ub
            xk_c = xk - inv(W) * D_c' * inv(D_c*inv(W)*D_c')*(D_c*xk-ub);
        end
        xk = xk_c;
        
    end
    
    % endpoint constraints
    if nargin > 14
        C1 = make_C(K,tdat(1),d,J);
        D_c = [C1 0 0];
        d_c = varargin{5};
        if nargin > 15
            C2 = make_C(K,tdat(end),d,J);
            D_c = [C1 0 0; C2 0 0];
            d_c = [varargin{5};varargin{6}];
        end
        W = inv(Pk);
        if nargin > 16 && varargin{7} == 1
            W = 1;
        end
        xk_c = xk - inv(W) * D_c' * inv(D_c*inv(W)*D_c')*(D_c*xk-d_c);
        xk = xk_c;
    end
    

    x1_n(:,k) = xk;
    P1_n(:,:,k) = Pk;
    
    hk = xk(1:J);
    
    zfit(k) = xk(J+1) + xk(J+2) * ytild(Ck*hk);
    phase1_n(k) = Ck*hk;
    
    xprev = xk;
    Pprev = Pk;
end

%hn_func = @(t) fnval(spmak(K,hk(1:J)'),t);
phase_func = @(t) fnval(spmak(K,hk'),t);
freq_func = @(t) fnval( fnder( spmak(K,hk')),t);

end












function [phi_func,spl,distort_vals] = ph_lin_distort(w,J,d,t,ph_noise,seed)
% create phase function that is near-linear with phase distortions
% (quasi-periodic)
% phi(t) = wt + f(t), f(0)<=a<=f(t)<=b
% phi'(t) = w + f'(t) >= 0 --> f'(t) >= -w
% Model distortions as B-spline function with derivative >= -w

Ktemp = -d:J;
K = (t(end)+.0001) * Ktemp / Ktemp(end-d);

% Randomly select B-spline coefficients based on dist_amp
rng(seed)
xcoefs = randn(1,J) * ph_noise;

% check if minimum derivative less than -w
spl = spmak(K,xcoefs);
splderv = fnder(spl);
minderv = min(fnval(splderv,t));
if minderv < -w
    error('set ph_noise smaller or w higher')
end

distort_vals = @(s) fnval(spl,s);
phi_func = @(s) w*s + fnval(spl,s);

end

function [x1_n,P1_n] = raw_phase_estimator(x0,P0,xdat,tdat,Q,R,ytild,dytild,varargin)

phase_ord = length(x0)-3;
if phase_ord > 2
    error('assume phase linear or quadratic (no higher orders for now)')
end

x1_n = zeros(4,length(tdat));
P1_n = [];
% raw_phase = zeros(1,length(tdat));
% freq = zeros(1,length(tdat));
dt = mean(diff(tdat));
if nargin > 9 && ~isempty(varargin{2})
    dt = varargin{2};
end
xprev = x0;
Pprev = P0;

ndat = length(xdat);
if nargin > 10
    ndat = varargin{3};
end

for k = 1:ndat
    zk = xdat(k);
    Ak = diag(ones(1,length(x0))); 
    if phase_ord == 1
        Ak(1,2) = dt;
    elseif phase_ord == 2
        Ak(1,2) = dt; Ak(1,3) = .5*dt^2;
        Ak(2,3) = dt;
    end
    xk_pred = Ak * xprev;
    
    Ck = [];
    innovk = [];
    if phase_ord == 1
        Ck = [xk_pred(3)*dytild(xk_pred(1)) 0 ytild(xk_pred(1)) 1];
        innovk = zk - (xk_pred(3)*ytild(xk_pred(1))+xk_pred(4));
    elseif phase_ord == 2
        Ck = [xk_pred(4)*dytild(xk_pred(1)) 0 0 ytild(xk_pred(1)) 1];
        innovk = zk - (xk_pred(4)*ytild(xk_pred(1))+xk_pred(5));
    end
    
    [xk,Pk] = kalman(xprev,Pprev,0,zk,Ak,0,Ck,Q,R,innovk);
    
    if nargin > 8
        if xk(2) < varargin{1}
            xk(2) = varargin{1};
           % xk(1) = xprev(1) + dt*xk(2);
%         elseif xk(1) < xprev(1)
%             xk(1) = xk_pred(1);
        end
    end
    
    x1_n(:,k) = xk;
    P1_n(:,:,k) = Pk;
    
    xprev = xk;
    Pprev = Pk;
    
end

end

function [phase_func,x1_n,P1_n,zfit] = phase_function_estimator(x0,P0,tdat,xdat,d,J,Q,R,ytild,dytild,varargin)

Ktemp = -d:J;
K = (tdat(end)+.001)*Ktemp/Ktemp(end-d);

xprev = x0;
Pprev = P0;

ndat = length(tdat);
if nargin > 11
    ndat = varargin{2};
end

x1_n = zeros(J+2,ndat);
P1_n = [];
zfit = zeros(1,ndat);

for k = 1:ndat
    [xk,Pk,~,~,~,~,~,~,~,Ck] = kalman_warp4(K,xprev,Pprev,R,tdat(k),xdat(k),Q,ytild,dytild);
    
    % increasing coefficients constraint
    if nargin > 10 && varargin{1} == 1
        for i = 2:J
            if xk(i) < 0
                xk(i) = 0;
            end
        end 
    end
    
    
    hk = [cumsum(xk(1:J)); xk(J+1:J+2)];
    
    x1_n(:,k) = xk;
    P1_n(:,:,k) = Pk;
    
    zfit(k) = xk(J+1) + xk(J+2) * w(Ck*hk);
        
    xprev = xk;
    Pprev = Pk;
end

hn_func = @(t) fnval(spmak(K,hk(1:J)'),t);
phase_func = hn_func;

end


function [phase_harm,a_harm,zfit] = get_phase_harmonics(x0,P0,Q,R,tdat,xdat,varargin)

r = (length(x0)-2)/2;
phase_harm = zeros(r,length(tdat));
zfit = zeros(1,length(tdat));

xprev = x0;
Pprev = P0;

for k = 1:length(tdat)
    
    dt = mean(diff(tdat));    
    
    Ak = diag(ones(1,2*r+2));
    for i = 1:r
        Ak(i,r+1) = i*dt;
    end
    xk_pred = Ak*xprev;
      
    Hk = zeros(1,2*r+2);
    for j = 1:length(Hk)
        if j <= r
            Hk(j) = xk_pred(j+r+1) * cos(xk_pred(j));
        elseif j > r+1 && j < 2*r+2
            Hk(j) = sin(xk_pred(j-(r+1)));
        elseif j == 2*r+2
            Hk(j) = 1;
        end
    end
    
    phik_pred = xk_pred(1:r);
    ak_pred = xk_pred(r+2:2*r+1);
    bk_pred = xk_pred(2*r+2);
    innovk = xdat(k) - (sum(ak_pred.*sin(phik_pred)+bk_pred));
    
    [xk,Pk] = kalman(xprev,Pprev,0,xdat(k),Ak,0,Hk,Q,R,innovk);
    
    if nargin > 6 && varargin{1} == 1
        if xk(r+1) < 0
            xk(r+1) = -xk(2);
            xk(1:r) = Ak(1:r,1:r+1)*xk(1:r+1);
        end
    end
    
    
    ak = xk(r+2:end-1);
    phik = xk(1:r);
    bk = xk(end);
    
    a_harm = ak;
    phase_harm(:,k) = phik;
    zfit(k) = sum(ak.*sin(phik)+bk);
    xprev = xk;
    Pprev = Pk;
end
end


function [f_qrs, df_qrs] = get_ecg_func(varargin)
fs2 = 256*10;

amp = [1.2 -5 30 -7.5 0.75];
if nargin > 0 && ~isempty(varargin{1})
    amp = varargin{1};
end

[s,ipeaks] = ecgsyn(fs2,30,0,60,0,0.5,fs2,[-70 -15 0 15 120],amp,[0.25 0.1 0.1 0.1 0.4]); %ai=[1.2 -5 30 -7.5 0.75] [0 0 30 -20 0]
rlocs = find(ipeaks==3);
qrs = s((rlocs(10):rlocs(11)-1));
qrs = circshift(qrs,round(length(qrs)/2));
t = (0:(length(qrs)-1))/fs2;
d = 3;
if nargin > 1
    d = varargin{2};
end
J = 100;
if nargin > 2
    J = varargin{3};
end
    
Ktemp = -d:J; K = Ktemp / Ktemp(end-d);
C = make_C(K,t,d,J);
xspl_coefs = inv(C'*C)*C'*qrs;
D_c = make_C(K,t(end),d,J) - make_C(K,t(1),d,J);
d_c = 0;
xspl_coefs_c = xspl_coefs - D_c' * inv(D_c*D_c') * (D_c*xspl_coefs - d_c);

spl = spmak(K,xspl_coefs_c');
spl_der = fnder(spl);

f_qrs = @(t) fnval(spl,mod(t,1)); %@(t) fnval(qrs_sp,mod(t-.5,1));
df_qrs = @(t) fnval(spl_der,mod(t,1));%@(t) fnval(qrs_sp_der,mod(t-.5,1));

end