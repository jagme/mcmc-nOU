function f = intensityFun(t, t0, period, delta, eta)
%INTENSITYFUNCTION Time-dependent jump intensity function of Geman and
%Roncoroni (2006)
%
% f = intensityFun(t,t0,period,delta,eta)
%
%         |                                 |
f = eta*( 2./ (1 + abs(sin(2*pi*(t-t0) / period)) ) -1).^delta;

%f = d + theta*cos(4*pi*(t-t0)/period);

%f = d*exp(theta*cos(4*pi*(t-t0)/period));