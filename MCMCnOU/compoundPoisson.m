function [X, T] = compoundPoisson(Tmax,t0,period,delta,eta,jParameters,jDist)
%COMPOUNDPOISSON Simulate paths of inhomogeneous Poisson process 
%with intensity function given in file intensityFun.m
%
% [X, T] = compoundPoisson(Tmax,t0,period,delta,eta,jParameters,jDist)
%
% Input arguments:
%
%   Tmax -  Maximum simulation time, ie. paths live on [0, Tmax]
%
%   t0, period, delta, eta are arguments to intensityFun.m
%
%   jParameters - parameters of the jump size distribution
%
%   jDist - jump size distribution as specified in drawJumpSize.m
%
% Output arguments: 
%
%   T - arrival times of the Poisson process path
%   X - corresponding jump sizes

lambda = eta + 0.1;     % fixed intensity bound, must be larger than lambda(t)
T = 0; X = 0;
t= exprnd(1/lambda);    % first time
i = 2; N = 0;           % number of jumps
% get jump times
while t <= Tmax
    U = rand;    %
    lambda_t = intensityFun(t,t0,period,delta,eta);
    if U <= lambda_t/lambda % accept
        N = N + 1;          % count jumps
        T(i) = t;
        X(i) = drawJumpSize(jParameters,jDist);
        i = i+1;        
    end
    t = t + exprnd(1/lambda);   % generate possible new jump time
end

