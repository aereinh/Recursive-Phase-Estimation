function [xk,Pk,hk,Fk_jac,Kk,innovk,resk,xpredk,Ppredk,Ck] = kalman_warp4(K,xprev,Pprev,Rk,sk,zk,Qk,w,dw,varargin)
% EKF recursion call to unwarp given template (w) to data (sk,zk)
% Assumes amplitude variation is linear (shifted and scaled)
% and phase variation is described by unconstrained B-Spline
%
%
% INPUTS: K -- knot vector (J+d+1)
%         xprev -- Previous (unconstrained) state estimate
%         Pprev -- Previous estimate covariance matrix
%         Rk -- Measurement noise covariance matrix
%         sk -- Current signal time
%         zk -- Current signal value
%         Qk -- State noise covariance matrix
%         w -- Template function (handle)
%         dw -- Template function derivative (handle)
%         [Optional inputs]
%         constraint -- 'unconstrained','exp','xsq', or 'logistic'
%         constraint_params -- [L,k,x0] for logistic
%         
%
% OUTPUTS: xk -- Posterior state estimate
%          Pk -- Estimate covariance matrix
%          hk -- Constrained state parameters 
%                 (1 spline function offset + J-1 spline coefficients + 2 amp params)
%          Fk_jac -- Jacobian for state-measurement transition
%          Kk -- Kalman gain
%          innovk -- Innovation (prediction vs. actual)
%          residk -- Residual (current estimate vs. actual)
%          xpredk -- Prediction for xk
%          Ppredk -- Prediction covariance
%          Ck -- B-spline evaluator matrix


% check inputs
if ~isa(w,'function_handle') || ~isa(dw,'function_handle')
    error('Input function handles for template function (w) and derivative (dw)')
end

% Obtain spline order (d) and number of coefficients (J)
J = length(xprev)-2;
nK = length(K);
d = nK-J-1;

% Define model components
Tild = diag(ones(1,J+2));
for i = 1:J
    Tild(i,1:i) = 1;
end
u = 0; % no input
B = 0; % no input
Ak = diag(ones(1,J+2)); % state transition matrix

% Use current timepoint (sk) to define new parameters
muk1 = find(K<=sk,1,'last'); % knot interval
% muk2 = find(K>=sk,1,'first');

muk = muk1;
% if ~isempty(muk2) && muk2 >= d+1
%     muk = min(muk1,muk2);
% end

% State transition (predict based on previous knots)
if nargin > 9 
    x0 = varargin{1};
    if ~isempty(x0)
        if xprev(muk) == x0(muk)
            %Ak(muk,2:muk-1) = 1/(muk-2);
            %Ak(muk,muk) = 0;
            Ak(muk,muk-1) = .9;
            Ak(muk,muk) = .1;
        end
    end
end
            
xk_p = Ak*xprev+B*u; % prediction (use to get innovation)
ampk_p = xk_p(J+1:J+2);
phasek_p = xk_p(1:J);

constraint = 'unconstrained';

if nargin > 10
    constraint = varargin{2};
    types = ["unconstrained","exp","xsq","logistic"];
    if sum(types==constraint) == 0
        error('constraint types: unconstrained, exp, xsq, logistic')
    end
end

%ampk_p(2) = exp(ampk_p(2));

if strcmp(constraint,'exp')
    ampk_p(2) = exp(ampk_p(2));
    phasek_p(2:J) = exp(phasek_p(2:J));
elseif strcmp(constraint,'xsq')
    ampk_p(2) = ampk_p(2)^2;
    phasek_p(2:J) = (phasek_p(2:J)).^2;
elseif strcmp(constraint,'logistic')
    L = 1; k = 1; c = 0;
    if nargin > 11
        params = varargin{3};
        L = params(1);
        k = params(2);
        c = params(3);
    end
    %ampk_p(2) = L/(1+exp(-k*ampk_p(2)))+c;
    phasek_p(2:J) = L/(1+exp(-k.*phasek_p(2:J)))+c;
end
    

hk_p = Tild*[phasek_p;ampk_p];


% Find Ck (b-spline evaluator)
bspline_vector = 1;
for j = d:-1:1
    B_term = make_B(muk,j,sk,K);
    bspline_vector = B_term*bspline_vector;
end
Ck = [zeros(1,(muk-(d+1))),bspline_vector,zeros(1,(J-muk))];

% innovation
innovk = zk - (ampk_p(1)+ampk_p(2)*w(Ck*hk_p(1:J)));

% Calculate Jacobian of state->measurement function
Dk = zeros(1,J);
for j = 1:J
    Dk(j) = sum(Ck(j:end));
end

Fk_phase = ampk_p(2)*dw(Ck*hk_p(1:J))*Dk;
Fk_amp = [1 w(Ck*hk_p(1:J))];
%Fk_amp(2) = Fk_amp(2)*exp(xk_p(J+2));

if strcmp(constraint,'exp')
    Fk_phase(2:J) = Fk_phase(2:J).*exp(xk_p(2:J))';
    Fk_amp(2) = Fk_amp(2)*exp(xk_p(J+2));
elseif strcmp(constraint,'xsq')
    Fk_phase(2:J) = Fk_phase(2:J).*xk_p(2:J)'*2;
    Fk_amp(2) = Fk_amp(2)*2*xk_p(J+2);
elseif strcmp(constraint,'logistic')
    dLogist_phase = (L*k*exp(-k*xk_p(2:J)'))/(1+exp(-k*xk_p(2:J)')).^2;
    %dLogist_amp = (L*k*exp(-k*xk_p(J+2)))/(1+exp(-k*xk_p(J+2)))^2;
    Fk_phase(2:J) = Fk_phase(2:J).*dLogist_phase;
    %Fk_amp(2) = Fk_amp(2)*dLogist_amp;
end
   
   
Fk_jac = [Fk_phase Fk_amp];


% Run EKF
[xk,Pk,xpredk,Ppredk,Kk] = kalman(xprev,Pprev,u,zk,Ak,B,Fk_jac,Qk,Rk,innovk);

phasek = xk(1:J);
ampk = xk(J+1:J+2);
%ampk(2) = exp(ampk(2));

if strcmp(constraint,'exp')
    ampk(2) = exp(ampk(2));
    phasek(2:J) = exp(phasek(2:J));
elseif strcmp(constraint,'xsq')
    ampk(2) = ampk(2)^2;
    phasek(2:J) = (phasek(2:J)).^2;
elseif strcmp(constraint,'logistic')
    %ampk(2) = L/(1+exp(-k*ampk(2)))+c;
    phasek(2:J) = L/(1+exp(-k.*phasek(2:J)))+c;
end

hk = Tild*[phasek;ampk];
resk = zk - (ampk(1)+ampk(2)*w(Ck*hk(1:J)));

end