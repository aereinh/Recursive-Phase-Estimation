function [Kt,xt_post,Pt_post,At,Qt,Ct,Rt] = kalman_bspline(Kt_prev,xt_prev,Pt_prev,Rt,st,yt,Ktbar,xbar,pbar,qbar)

J = length(xt_prev);
K = length(Kt_prev);
d = K-J-1;
I = J-d;
sigma = 0;

% Check if knot shift is necessary
if st >= Kt_prev(J+1)
    if st < Kt_prev(1)
        sigma = -(d+1);
    else
        % find sigma s.t. st in [Kt_prev(d+1+sigma),Kt_prev(d+2+sigma))
        index = find(Kt_prev<=st,1,'last');
        sigma = index-d-1;
    end
end

% Set At based on Eq (16)
At = zeros(J,J);
for g = 1:J
    for h = 1:J
        if h == (g + sigma)
            At(g,h) = 1;
        end
    end
end

Qt = qbar*diag(ones(1,J));

if sigma >= 0
    % Set Kt as in Eq (15)
    Kt = [Kt_prev((sigma+1):K), Ktbar((K-sigma+1):K)];
    ut = [zeros(1,J-sigma),xbar*ones(1,sigma)]';
    for m = J-sigma+1:J
        Qt(m,m) = pbar;
    end
else
    % Set Kt as in Eq (19)
    Kt = [Ktbar(1:(-sigma)), Kt_prev(1:(K+sigma))];
    ut = [xbar*ones(1,-sigma),zeros(1,J+sigma)]';
    for m = 1:(-sigma)
        Qt(m,m) = pbar;
    end
end

% Set mu s.t. st in [Kt(mu),Kt(mu+1))
mu = find(Kt<=st,1,'last');
if st >= Kt(mu+1)
    error('Could not find proper mu')
end

Vt = length(yt);

% Set Ct in R(VtxJ) as in Eq (20)
Ct = zeros(Vt,J);
for v = 1:Vt
    r = v-1; % Assumes measurements in vth row refer to (v-1)th derivative
    bspline_vector = zeros(1,(d+1));
    if r <= d
        bspline_vector = 1;
        for j = 0:(r-1)
            B_prime_term = make_B_prime(mu,d-j,Kt);
            bspline_vector = B_prime_term*bspline_vector;
        end
        for k = (d-r):-1:1
            B_term = make_B(mu,k,st,Kt);
            bspline_vector = B_term*bspline_vector;
        end
        bspline_vector = factorial(d)/factorial(d-r) * bspline_vector;
    end
    c = [zeros(1,(mu-(d+1))),bspline_vector,zeros(1,(J-mu))];
    Ct(v,:) = c;
end


ident = diag(ones(1,J));

[xt_post,Pt_post] = kalman(xt_prev,Pt_prev,ut,yt,At,ident,Ct,Qt,Rt);


end