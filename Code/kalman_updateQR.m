function [Qk,Rk] = kalman_updateQR(Qkprev,Rkprev,innk,resk,Ck,P_predk,Kk,alphaq,alphar)
%Updates state and measurement noise covariance matrices
% INPUTS: 
%   Qkprev
%   Rkprev
%   innk
%   resk
%   Ck
%   P_predk
%   Kk
%   alphaq
%   alphar
%
% OUTPUTS
%   Qk
%   Rk


Qk = alphaq*Qkprev + (1-alphaq)*(Kk*(innk*innk')*Kk');
Rk = alphar*Rkprev + (1-alphar)*(resk*resk'+Ck*P_predk*Ck');

end

