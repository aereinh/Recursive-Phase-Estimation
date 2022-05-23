function B_prime = make_B_prime(mu,delta,K)

B_prime = zeros(delta,(delta+1));

for i = 1:delta
    B_prime(i,i) = -1/(K(mu+i)-K(mu+i-delta));
    B_prime(i,i+1) = 1/(K(mu+i)-K(mu+i-delta));
end

end