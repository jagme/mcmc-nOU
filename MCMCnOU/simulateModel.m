function [E,T,L] = simulateModel(Nt, mu, lambda0, sigma, lambda1, jRate, beta, negJ)
%SIMULATEMODEL Simulate paths of the 2-OU models, X = Y0 + Y1
%
% [E,T,L] = simulateModel(Nt, mu, lambda0, sigma, lambda1, jRate, beta, negJ)
%
% dY0 = 1/lambda0*(mu-Y0)dt + sigma*dWt, Wt is a Wiener process
%
% dY1 = -1/lambda1*Y1*dt + dL1, L1 is a homogeneous or inhomgeenous 
% compound Poisson process
%
% Input arguments:
%
%   Nt - Number of timesteps, observation period is [0,Nt] with grid size 1
%
%   mu - Level of mean reversion of Y0
%
%   lambda0 - Time to mean reversion (inverse of speed of mean reversion)
%
%   sigma - Volatility
%
%   lambda1 - Time to mean reversion of Y1
%
%   jRate - if L1 is time-homogeneous, jRate is a scalar representing the 
%       jump intensity rate of L1. If L1 is time-inhomogeneous jRate is 
%       structure array containing the parameters for the intensity
%       function given in intensityFun.m, eg. jRate.t0 = 130, 
%       jRate.period = 260, jRate.delta = 0.15, jRate.eta = 0.3
%
%   beta - Mean jump size of L1
%
%   negJ  - scalar indicating whether Y1 is positive (negJ = 1) 
%               or negative (negJ = -1)
%
% Output arguments:
%
%   E - Discretised sample path at grid times 0,1,...,Nt, E1 below is
%       the exact solution
%   T - True jump times
%   L - True jump sizes

% choose jump size distribution
%   1 for Pareto, 2 Exponential, 3 Gamma
jumpDistribution = 2;

% simulate Poisson process

if ~isstruct(jRate)   % homogeneous
    
    [L,jumpTimes] = homogeneousCPoisson(Nt, jRate, beta, jumpDistribution);
    
else                  % inhomogeneous
    
    [L,jumpTimes] = compoundPoisson(Nt, jRate.t0, jRate.period, ...
        jRate.delta, jRate.eta, beta, jumpDistribution);
    
end
% simulation times, combine fix time grid points and jump times
[T,IX] = sort([1:Nt jumpTimes]);    % IX is the permutation vector
L = [zeros(1,Nt) L]; L = L(IX);     % sort L over the 'combined' grid

lenT = length(T); % time

% initialise variables
E1 = zeros(1,lenT); 
Y0 = zeros(1,lenT); 
Y1 = zeros(1,lenT);
E1(1) = mu; Y0(1) = E1(1); % initial values at time 0
%dt = 1;

% simulate path
for i = 1:lenT-1
    dt = T(i+1) - T(i); % dt can change from step to step
    temp = exp(-1/lambda0*dt);
    Y0(i+1) = mu*(1-temp) +  temp*Y0(i) + sigma*sqrt((1-temp^2)/2*lambda0)*randn;    
    Y1(i+1) = exp(-1/lambda1*dt)*Y1(i) + L(i+1);
    E1(i+1) = Y0(i+1) + negJ*Y1(i+1);
end

% E1 contains exact simulation of the process
% the simulation times are contained in T
% E below contains the observations at times 0, 1, 2, ..., Nt

E = 0; j = 0;
for i=1:length(T)
    if T(i)==j
        E(j+1) = E1(i);
        j = j+1;
    end
end

