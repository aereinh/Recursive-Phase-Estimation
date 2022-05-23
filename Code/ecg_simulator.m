% Creates synthetic ecg signal

% Set parameters
thetai = [-1/3 * pi, -1/12 * pi, 0, 1/12 * pi, 1/2 * pi];
ai = [1.2, -5, 30, -7.5, .75];
bi = [.25, .1, .1, .1, .4];
w = 1;

alpha = @(x,y) 1-sqrt(x.^2+y.^2);
theta = @(x,y) atan2(y,x);
xdot = @(t,W) alpha(W(1),W(2))*W(1)-w*W(2);
ydot = @(t,W) alpha(W(1),W(2))*w*W(1);
zdot = @(t,W) - sum(ai.*mod(theta(W(1),W(2))-thetai,2*pi).*exp(mod(theta(W(1),W(2))-thetai,2*pi).^2./(2*bi)) - W(3));
Wdot = @(t,W) [xdot(t,W); ydot(t,W); zdot(t,W)];
W0 = [0;0;1];

SOL = ode45(Wdot,[-1 1],W0);

tvals = linspace(-1,1,1000);
wvals = deval(SOL,tvals);


figure(1)
clf
plot(tvals,wvals(3,:))