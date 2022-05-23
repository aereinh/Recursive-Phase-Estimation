%% Load in signal

load('ECGSignal')
ecg = ECGSignal(:,1);
L = length(ecg);
t = 0:100/L:100-100/L;
fs = 1/mean(diff(t));

L_samp = 5000;
ecg_samp = ecg;
t_samp = t(1:L_samp);
ecg_filt = lowpass(ecg_samp,10,fs);

[pks,locs] = findpeaks(ecg_filt,'MinPeakHeight',.2);
max_interval = max(diff(locs));

ecg_frames = NaN(max_interval,length(pks)-1);

for i = 1:length(pks)-1
    frame_len = locs(i+1)-locs(i)+1;
    ecg_frames(1:frame_len,i) = ecg_filt((locs(i)-locs(1)+1):(locs(i+1)-locs(1)+1));
end

ecg_frames(end,:) = NaN;

%%
sig = ecg_filt'-mean(ecg_filt);
st = t;

% Use spectrum to find fundamental frequency

spect = fft(sig);
L = length(spect);
spect = spect(1:L/2+1);
spect(2:end) = 2*spect(2:end);
fz = linspace(0,fs/2,length(spect));

pk_index = find(abs(spect) == max(abs(spect))); % only works for relatively clean data
ffreq = fz(abs(spect) == abs(spect(pk_index)));
phase_offset = angle(spect(pk_index));

sin_func = @(t) sin(2*pi*ffreq*t+phase_offset);
sin_vals = sin_func(st);

figure(1)
clf
hold on
plot(st,sig)
plot(st,sin_vals)


%% Use KF to estimate amplitude, frequency, and phase offset of sine wave in noise 
% (no time variance in parameter value)
rng(1234)

% atrue = 1;
% wtrue = pi;
% thetatrue = -pi/2;

w0 = @(t,a,w,theta) a*cos(w*t+theta);

% st = linspace(0,10,2000);
% zt = f(st,atrue,wtrue,thetatrue)+0.05*randn(size(st));
st = t;
zt = sig;


x0 = [1;1;0]; % Initial guess
P0 = 10^(2)*diag([1,1,1]); % Initial state (co)variance
u = 0; % Assume no input signal
B = 1; % Input matrix
A = diag([1,1,1]); % State transition matrix
Q = 10^(-12)*diag([1,1,1]); % State covariance matrix
R = 1;

% Run extended KF for first data point
s1 = st(1);
z1 = zt(1);
xprev = x0;
Pprev = P0;
x1_pred = A*xprev+B*u;
a1p = x1_pred(1);
w1p = x1_pred(2);
t1p = x1_pred(3);
C1 = [sin(w1p*s1+t1p), a1p*s1*cos(w1p*s1+t1p), a1p*cos(w1p*s1+t1p)]; % Jacobian
res = z1-a1p*sin(w1p*s1+t1p);
[x1,P1] = kalman(xprev,Pprev,u,z1,A,B,C1,Q,R,res);

% Set resultant values
x1_n = x1;
P1_n = P1;

for k = 2:length(st)
    sk = st(k);
    zk = zt(k);
    xprev = x1_n(:,k-1);
    Pprev = P1_n(:,:,k-1);
    xk_pred = A*xprev+B*u;
    akp = xk_pred(1);
    wkp = xk_pred(2);
    tkp = xk_pred(3);
    Ck = [sin(wkp*sk+tkp), akp*sk*cos(wkp*sk+tkp), akp*cos(wkp*sk+tkp)];
    res = zk-akp*sin(wkp*sk+tkp);
    [xk,Pk] = kalman(xprev,Pprev,u,zk,A,B,Ck,Q,R,res);
    x1_n(:,k) = xk;
    P1_n(:,:,k) = Pk;
end

x1_n(:,end)./[1;(2*pi);1]
P1_n(:,:,end)


%% 

aest = x1_n(1,end);
west = x1_n(2,end);
test = x1_n(3,end);

f1vals = w0(st,aest,west,test);


figure(1)
clf
hold on
plot(st,zt)
plot(st,25*f1vals)

%% Time-varying phase (B-spline)
rng(1234)

% variables to adjust model
Jtemp = 20;
P0temp = 10^(4);
Rtemp = 25;
noise = .02;
lbtemp = 0.1; % between 0 and 1
uptemp = inf; % above 1
ntemp = 1; % between 0 and 1
constrained = 0; % for computing constrained estimates
constrained2 = 0; % for feeding constrained estimate into EKF
constrained_plot = 0; % for plotting constrained estimates
typeW = 'ls'; % ls or ml
dK_adj1 = 0; % 0 or 1
dK_adj2 = 0; % 0 or 1


% Simulated data
htrue = @(t) t.^5;
wwarped = @(t) sin(pi*htrue(t));
tdat = linspace(0,1,100);
tdat(end) = [];
xdat = wwarped(tdat)+noise*randn(1,length(tdat));
ndat = length(tdat);

% (Un)warp a sine template to the data directly (i.e. no state function estimates)

d = 3; % degree
J = Jtemp; % number of coefficients
nK = J+d+1; % number of knots
Ktemp = ((1:nK)-4); % knots (uniform spacing)
K = Ktemp / Ktemp(end-d);
dK = mean(diff(K));
tau = 1/Ktemp(end-d); % spacing in hk coefficients to get [K(d+1),K(J+1)]->[0,1] mapping

w0 = @(t) sin(pi/K(end-d)*t); % template function
dw0 = @(t) pi/K(end-d) * cos(pi/K(end-d)*t);

% New state variable hk_tild is consecutive differences of hk elements
% This can now be constrained to be positive

% Set initial values

g0 = ones(J-1,1)*tau;
P0 = P0temp * diag(ones(1,J-1));

Tild = zeros(J,J-1);
for i = 2:J
    Tild(i,1:(i-1)) = ones(1,i-1);
end
h0 = Tild*g0;
h0_func = @(t) fnval(spmak(K,h0'),t) - fnval(spmak(K,h0'),K(4)); % initial warping function

mu0 = d+1;
Pprev = P0;
gprev = g0;
hprev = h0;


% Set other KF parameters
u = 0; % Assume no input signal
B = 1; % Input matrix
A = 1;% A = diag(ones(1,J)); % State transition matrix
Q = 10^(-12)*diag(ones(1,J-1)); % State covariance matrix
R = Rtemp;
A_constr = diag(-ones(1,J-1)); A_constr(1:d+1,1:d+1) = 0; % constraint matrix
b_constr = -lbtemp*tau*ones(J-1,1); b_constr(1:d+1) = 0; % constraint value
% lowerbound = lbtemp*tau*ones(J-1,1);
% upperbound = uptemp*tau*ones(J-1,1);
ginit = ones(J-1,1); % constraint initial search value

% Set reference matrix measurement (so that h(t) starts at 0 at K(d+1))
r = 0;
bspline_vector = 1;
for k = (d-r):-1:1
    B_term = make_B(mu0,k,K(d+1),K);
    if dK_adj1 == 1
        B_term = make_B(mu0,k,K(d+1)-dK,K);
    end
    bspline_vector = B_term*bspline_vector;
end
Cref = [zeros(1,(mu0-(d+1))),bspline_vector,zeros(1,(J-mu0))];

% Set matrices of estimates and covariances
g1_n_unconstr = [];
g1_n_constr = [];
h1_n_unconstr = [];
h1_n_constr = [];
P1_n = [];

% Predict and update using new data
for k = 1:ndat
    sk = tdat(k);
    zk = xdat(k);
    muk = find(K<=sk,1,'last');
    gk_p = A*gprev+B*u;
    hk_p = Tild*gk_p;
    
    % Matrix for evaluating h(s)
    r = 0;
    bspline_vector = 1;
    for j = (d-r):-1:1
        B_term = make_B(muk,j,sk,K);
        if dK_adj1 == 1
            B_term = make_B(muk,j,sk-dK,K);
        end
        bspline_vector = B_term*bspline_vector;
    end
    Ck = [zeros(1,(muk-(d+1))),bspline_vector,zeros(1,(J-muk))];
    Ck = Ck-Cref;

    % Measurement matrix (derivative of w0(Ck*Tild*g)) evaluated at g1_p
    Cksum = zeros(1,J-1);
    Cksum(1:(muk-4)) = sum(Ck((muk-3):muk));
    Cksum(muk-3) = sum(Ck((muk-2):muk));
    Cksum(muk-2) = sum(Ck((muk-1):muk));
    Cksum(muk-1) = Ck(muk);
    
    Dk = Cksum * dw0(Ck*hk_p);
    
    % Dk = Ck * dw0(Ck*hk_p);
    
    resk = zk - w0(Ck*hk_p);
    
    [gk_unconstr,Pk] = kalman(gprev,Pprev,u,zk,A,B,Dk,Q,R,resk);
    hk_unconstr = Tild*gk_unconstr;
    
    g1_n_unconstr(:,k) = gk_unconstr;
    h1_n_unconstr(:,k) = hk_unconstr;
    P1_n(:,:,k) = Pk;
    
    % Solve inequality constraint problem if necessary
%     if any(gk_unconstr<lbtemp*tau) || any(gk_unconstr>uptemp*tau)

    gk_constr = gk_unconstr;
    if constrained == 1
        if any(A_constr*gk_unconstr > b_constr)
            if strcmp(typeW, 'ls')
                Wk = 1;
            elseif strcmp(typeW,'ml')
                Wk = Pk;
            end
            f_constr_k = @(g) (g-gk_unconstr)'*Wk*(g-gk_unconstr);
            gk_constr = fmincon(f_constr_k,ginit,A_constr,b_constr);
    %         gk_constr = fmincon(f_constr_k,ginit,[],[],[],[],lowerbound,upperbound);
        end
    end
    
    hk_constr = Tild*gk_constr;
    
    g1_n_constr(:,k) = gk_constr;
    h1_n_constr(:,k) = hk_constr;
    
    
    if constrained2 == 1
        gprev = gk_constr;
    else
        gprev = gk_unconstr;
    end
    Pprev = Pk;
    
    
end

% check_inc = zeros(1,size(h1_n,2));
% for i = 1:size(h1_n,2)
%     check_inc(i) = issorted(h1_n(:,i));
% end

traces = zeros(1,ndat);
for i = 1:ndat
    traces(i) = trace(P1_n(:,:,i));
end

% Plot results

n = round(ndat*ntemp);

hn = h1_n_unconstr(:,n);
if constrained_plot == 1
    hn = h1_n_constr(:,n);
end

hn_func = @(t) fnval(spmak(K,hn'),t) - fnval(spmak(K,hn'),K(4));
if dK_adj2 == 1
    hn_func = @(t) fnval(spmak(K,hn'),t-dK) - fnval(spmak(K,hn'),K(4)-dK);
end

tvals = linspace(K(d+1),K(end-d),100);
h0vals = h0_func(tvals);
hnvals = hn_func(tvals);
wvals = w0(tvals);
wwarpedvals = wwarped(tvals);
% xvals = g(tvals);

figure(1)
clf
subplot(411)
plot(tvals,wvals)
title('template')
subplot(412)
hold on
plot(tvals,htrue(tvals))
plot(tvals,h0vals)
plot(tvals,hnvals)
title('warping functions')
subplot(413)
hold on
plot(tvals,w0(h0vals),'b')
plot(tvals,w0(hnvals),'r')
plot(tdat,xdat,'black.')
plot(tdat(n),xdat(n),'r*')
plot(tvals,wwarpedvals)
title('warped template')
subplot(414)
plot(tdat,traces)
title('trace')


%% **** EKF for recursive estimation of warping functions 
% Try putting nonnegativity constraints using another state

%rng(1235)

% Simulated data
anoise = .5;
beta = -2; T = 1/.5;
%htemp = @(t) exp(3*t);
trueh = @(t) T*(exp(beta*t)-1)/(exp(beta*T)-1);
wwarped = @(t) sin(2*pi*trueh(t));
tdat = linspace(0,1,2000);
tdat(end) = [];
xdat = wwarped(tdat)+anoise*randn(1,length(tdat));
ndat = length(tdat);

% (Un)warp a sine template to the data directly (i.e. no state function estimates)
d = 3; % degree
J = 10; % number of coefficients
nK = J+d+1; % number of knots
Ktemp = ((1:nK)-4); % knots (uniform spacing)
K = Ktemp / Ktemp(end-d);
tau = mode(diff(K));

% Define template and template derivative
w0 = @(t) sin(2*pi*t); % template function
dw0 = @(t) 2*pi * cos(2*pi*t);

% Set initial values
x0 = zeros(J-1,1);
g0 = exp(x0) * tau ;
Tild = zeros(J,J-1);
for i = 2:J
    Tild(i,1:(i-1)) = ones(1,i-1);
end
h0 = Tild*g0;
h0_func = @(t) fnval(spmak(K,h0'),t) - fnval(spmak(K,h0'),K(4)); % initial warping function
pbar = 10^2;
qbar = 10^(-12);
Q = diag(qbar*ones(J-1));
R = 50;
P0 = pbar*diag(ones(1,J-1));
Pprev = P0;
xprev = x0;


% Initialize matrices to store results
x1_n = [];
g1_n = [];
h1_n = [];
P1_n = [];
F1_n = [];
eigs1_n = [];

ntemp = ndat;

for k = 1:ntemp %ndat
    [xk,Pk,gk,hk,~,~,~,Fk,Kk] = kalman_warp(K,xprev,Pprev,R,tdat(k),xdat(k),Q,w0,dw0);
    eigs = eig(diag(ones(1,length(Kk)))-Kk*Fk);
    x1_n(:,k) = xk;
    g1_n(:,k) = gk;
    h1_n(:,k) = hk;
    P1_n(:,:,k) = Pk;
    F1_n(k,:) = Fk;
    eigs1_n(:,k) = eigs;
    xprev = xk;
    Pprev = Pk;
end

traces = zeros(1,ndat);
for i = 1:ntemp
    traces(i) = trace(P1_n(:,:,i));
end

% Plot results
n = ntemp;
hn = h1_n(:,n);
hn_func = @(t) fnval(spmak(K,hn'),t) - fnval(spmak(K,hn'),K(4));
% k = round(n/2);
hk = h1_n(:,k);
hk_func = @(t) fnval(spmak(K,hk'),t) - fnval(spmak(K,hk'),K(4));
tvals = linspace(K(d+1),K(end-d),100);
h0vals = h0_func(tvals);
hnvals = hn_func(tvals);
% hkvals = hk_func(tvals);
htruevals = trueh(tvals);
wvals = w0(tvals);
wwarpedvals = wwarped(tvals);
% xvals = g(tvals);

figure(1)
clf
subplot(411)
plot(tvals,wvals)
legend('w(t)')
title('template')
subplot(412)
hold on
plot(tvals,h0vals)
% plot(tvals,hkvals)
plot(tvals,hnvals,'r')
plot(tvals,htruevals)
legend('h0(t)','hn(t)','h(t)')
title('warping functions')
subplot(413)
hold on
%plot(tvals,w0(h0vals))
% plot(tvals,w0(hkvals))
plot(tvals,w0(hnvals),'r','LineWidth',5)
plot(tdat,xdat,'black.')
plot(tdat(n),xdat(n),'g*')
% plot(tdat(k),xdat(k),'r*')
% plot(tvals,wwarpedvals)
%legend('w(h0(t))','w(hn(t))')
title('warped template')
subplot(414)
plot(tdat,traces)
title('trace(Pk)')

%% Include phase offset (kalman_warp3)

% Simulated data
anoise = .1;
beta = -2; T = 1/.5;
%htemp = @(t) exp(3*t);
%trueh = @(t) T*(exp(beta*t)-1)/(exp(beta*T)-1)+.5;
trueh = @(t) 2*t;
trueamp0 = 0.1;
trueamp1 = 1.5;
wwarped = @(t) trueamp1*sin(2*pi*trueh(t)) + trueamp0;
tdat = linspace(0,1,2000);
tdat(end) = [];
xdat = wwarped(tdat)+anoise*randn(1,length(tdat));
ndat = length(tdat);

% (Un)warp a sine template to the data directly (i.e. no state function estimates)
d = 3; % degree
J = 50; % number of coefficients
nK = J+d+1; % number of knots
Ktemp = ((1:nK)-4); % knots (uniform spacing)
K = Ktemp / Ktemp(end-d);
tau = mode(diff(K));

% Define template and template derivative
w0 = @(t) sin(2*pi*t); % template function
dw0 = @(t) 2*pi * cos(2*pi*t);

% Set initial values
x0 = zeros(J+2,1);
g0 = x0;
g0(2:J) = exp(x0(2:J))*tau;
g0(J+2) = exp(x0(J+2));
h0 = g0;
h0(2:J) = cumsum(g0(2:J));
h0_func = @(t) h0(1) + fnval(spmak(K,[0 h0(2:J)']),t) - fnval(spmak(K,[0 h0(2:J)']),K(d+1)); % initial warping function

% Noise parameters
pph0 = 1e-10;
pph1 = 1e1;
pamp0 = 1e1;
pamp1 = 1e1;
P0 = diag([pph0, repmat(pph1,1,J-1), pamp0, pamp1]);
qbar = 10^(-12);
Q = diag(qbar*ones(1,J+2));
R = 20;

Pprev = P0;
xprev = x0;

% Initialize matrices to store results
x1_n = [];
h1_n = [];
P1_n = [];
F1_n = [];
eigs1_n = [];

ntemp = ndat;

for k = 1:ntemp %ndat
    [xk,Pk,hk,Fk,Kk,innovk,resk,xpredk,Ppredk] = kalman_warp3(K,xprev,Pprev,R,tdat(k),xdat(k),Q,w0,dw0,tau);  
    eigs = eig(diag(ones(1,length(Kk)))-Kk*Fk);
    x1_n(:,k) = xk;
    h1_n(:,k) = hk;
    P1_n(:,:,k) = Pk;
    F1_n(k,:) = Fk;
    eigs1_n(:,k) = eigs;
    xprev = xk;
    Pprev = Pk;
end

traces = zeros(1,ndat);
for i = 1:ntemp
    traces(i) = trace(P1_n(:,:,i));
end

% Plot results
n = ntemp;
hn = h1_n(:,n);
hn_func = @(t) hn(1) + fnval(spmak(K,[0 hn(2:J)']),t) - fnval(spmak(K,[0 hn(2:J)']),K(d+1));
% k = round(n/2);
tvals = linspace(K(d+1),K(end-d),100);
h0vals = h0_func(tvals);
hnvals = hn_func(tvals);
% hkvals = hk_func(tvals);
htruevals = trueh(tvals);
wvals = w0(tvals);
wwarpedvals = wwarped(tvals);
% xvals = g(tvals);

figure(1)
clf
subplot(411)
plot(tvals,wvals)
legend('w(t)')
title('template')
subplot(412)
hold on
plot(tvals,h0vals)
% plot(tvals,hkvals)
plot(tvals,hnvals,'r')
plot(tvals,htruevals)
legend('h0(t)','hn(t)','h(t)')
title('warping functions')
subplot(413)
hold on
%plot(tvals,w0(h0vals))
% plot(tvals,w0(hkvals))
plot(tvals,hn(J+1)+hn(J+2)*w0(hnvals),'r','LineWidth',5)
plot(tdat,xdat,'black.')
plot(tdat(n),xdat(n),'g*')
% plot(tdat(k),xdat(k),'r*')
% plot(tvals,wwarpedvals)
%legend('w(h0(t))','w(hn(t))')
title('warped template')
subplot(414)
plot(tdat,traces)
title('trace(Pk)')
