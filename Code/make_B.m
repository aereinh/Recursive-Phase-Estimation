function B = make_B(mu,delta,s,K)

B = zeros(delta,(delta+1));

for i = 1:delta
    B(i,i) = (K(mu+i)-s)/(K(mu+i)-K(mu+i-delta));
    B(i,i+1) = (s-K(mu+i-delta))/(K(mu+i)-K(mu+i-delta));
end

end