function [x_post,P_post,x_pred,P_pred,K] = kalman(xprev,Pprev,u,y,A,B,C,Q,R,varargin)

% INPUTS: 
%   xprev -- Previous state estimate
%   Pprev -- Previous state estimate covariance
%   u -- Input/control vector
%   y -- Current observed value
%   A -- State transition matrix (or Jacobian of state transition function)
%   B -- Input/control matrix
%   C -- Measurement matrix (or Jacobian of measurement function)
%   Q -- State noise covariance matrix
%   R -- Measurement noise covariance matrix
%   [optional...] 
%   innov -- For EKF, input innovation using nonlinear measurement function
%   xpred -- For EKF, input prediction using nonlinear state transition
%   function
%
% OUTPUTS:
%   x_post -- Posterior state estimate 
%   P_post -- Posterior state estimate covariance
%   x_pred -- Predicted state estimate
%   P_pred -- Predicted state estimate covariance
%   K -- Optimal Kalman gain matrix


% Time update
x_pred = A*xprev+B*u;
% For EKF, Jacobian (A) of state transition function and nonlinear
% prediction is inputted
if nargin > 10
    if ~isempty(varargin{2})
        x_pred = varargin{2};
    end
end

P_pred = A*Pprev*A'+Q;

% Measurement update
S = C*P_pred*C'+R;
K = P_pred*C'*inv(S);
dim = size(K*C,1);
ident = diag(ones(1,dim));

innov = (y-C*x_pred);
% For EKF, Jacobian (C) of measurement function and nonlinear innovation is inputted
if nargin > 9
    if ~isempty(varargin{1})
        innov = varargin{1};
    end
end

x_post = x_pred+K*innov;
P_post = (ident-K*C)*P_pred*(ident-K*C)'+K*R*K';


end







