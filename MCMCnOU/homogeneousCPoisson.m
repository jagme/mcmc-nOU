function [X,T,n] = homogeneousCPoisson(Tmax,lambda,parameters,jDist)
%HOMOGENOUSCPOISSON Simulate homogeneous compound Poisson process
%
% [X,T,n] = homogeneousCPoisson(Tmax,lambda,parameters,jDist)
%
% Input arguments:
%
%   Tmax  - maximum simulation time, ie, paths live on [0, Tmax]
%
%   lambda - jump intensity per unit time, over [0 Tmax]
%
%   parameters - parameters for jump size disribution 
%
%   jDist - jump size distribution as specified in drawJumpSize.m           
%
% Ouput arguments: 
% 
%   T  - arrival times of the Poisson process
%
%   X  - corresponding jump sizes
%
%   n  - total number of jumps
%

T(1) = 0; % initial time
X(1) = 0;

t = exprnd(1/lambda); % possibly next jump time
i = 2;
% get jump times
while t <= Tmax    
    T(i) = t;      
    X(i) = drawJumpSize(parameters,jDist);
    i = i + 1;
    t = t + exprnd(1/lambda); % time of next jump
end
n = i - 2; % number of jumps

% plot compound Poisson process
%stairs(T,cumsum(X))