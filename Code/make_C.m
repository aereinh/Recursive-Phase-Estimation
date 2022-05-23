function C = make_C(K,t,d,J,varargin)
% INPUTS: K, t, d, J
% Makes B-spline evaluator matrix for knots K, time points t, degree d, and
% J coefficients
    C = zeros(length(t),J);
    
    for k = 1:length(t)
        muk = find(K>t(k),1,'first')-1;
        if muk >= J+1
            muk = J;
        end
        
        bspline_vector = 1;
        for j = d:-1:1
            B_term = make_B(muk,j,t(k),K);
            bspline_vector = B_term*bspline_vector;
        end
        
        % make derivative evaluator
        if nargin > 4
            r = varargin{1};
            bspline_vector = zeros(1,(d+1));
            if r <= d
                bspline_vector = factorial(d) / factorial(d-r);
                for delta = d:-1:(d-r+1)
                    bspline_vector = make_B_prime(muk,delta,K) * bspline_vector;
                end
                for delta = (d-r):-1:1
                    bspline_vector = make_B(muk,delta,t(k),K) * bspline_vector;
                end        
            end
        end 
        
        C(k,:) = [zeros(1,(muk-(d+1))),bspline_vector,zeros(1,(J-muk))];   
        
    end
    
        
    %D(k,:) = dK*flip(cumsum(flip(C(k,:))));
    
    
end