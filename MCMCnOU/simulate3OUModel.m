function [E,LT1,LT2] = simulate3OUModel(Nt, mu, lambda0, sigma, lambda1,lambda2, jRate1, jRate2, beta1, beta2, negJ)
%SIMULATE3OUMODEL Simulate the 3-OU model X = Y0 + w1*Y1 + w2*Y2
%
% [E,LT1,LT2] = simulate3OUModel(Nt, mu, lambda0, sigma, lambda1,lambda2, jRate1, jRate2, beta1, beta2, negJ)
%
% dY0 = 1/lambda0*(mu-Y0)dt + sigma*dWt, Wt is a Wiener process
%
% dYi = -1/lambdai*Yi*dt + dLi, Li is a compound Poisson process with
% jump intensity parameters jRatei and Exp(betai) distributed jump sizes
%
% Input arguments:
%
%   Nt      - Number of timesteps, observation period is [0,Nt] with grid size 1
%
%   mu      - Level of mean reversion of Y0
%
%   lambda0 - Time to mean reversion (inverse of speed of mean reversion)
%   of Y0
%
%   sigma   - Volatility
%
%   lambda1 - Time to mean reversion of Y1
%
%   lambda2 - Time to mean reversion of Y2
%
%   jRatei  - if Li is time-homogeneous, jRatei is a scalar representing 
%       the jump intensity rate of Li. If Li is time-inhomogeneous jRatei is 
%       a structure array containing the parameters for the intensity
%       function given in intensityFun.m, eg. jRatei.t0 = 130, 
%       jRatei.period = 260, jRatei.delta = 0.15, jRatei.eta = 0.3
%
%   beta1   - Mean jump size of L1
%
%   beta2   - Mean jump size of L2
%
%   negJ  - 2-by-1 vector indicating whether Yi is positive (negJ(i) = 1) 
%               or negative (negJ(i) = -1)
%
% Output arguments:
%
%   E       - Discretised sample path at grid times 0,1,...,Nt, E1 below is
%             the exact solution
%   LT1     - True jump process L1
%   LT2     - True jump process L2
%

% jump size distribution: 1 Pareto, 2 exponential, 3 Gamma
jumpDistribution = 2;
% generate jump times and jump sizes
if ~isstruct(jRate1)   % homogeneous    
    [L1,jumpTimes1] = homogeneousCPoisson(Nt,jRate1,beta1,jumpDistribution);    
else                  % inhomogeneous
    [L1,jumpTimes1]   = compoundPoisson(Nt, jRate1.t0, jRate1.period, ...
        jRate1.delta, jRate1.eta, beta1, jumpDistribution);  
end

if ~isstruct(jRate2)   % homogeneous
    [L2,jumpTimes2] = homogeneousCPoisson(Nt,jRate2,beta2,jumpDistribution);
else
    [L2,jumpTimes2]   = compoundPoisson(Nt, jRate2.t0, jRate2.period, ...
        jRate2.delta, jRate2.eta, beta2, jumpDistribution);  
end
% remove initial time 0
L1 = L1(2:end); jumpTimes1 = jumpTimes1(2:end); L2 = L2(2:end); jumpTimes2 = jumpTimes2(2:end); 
% simulation times, combine fix intervals and jump times
[T,IX] = sort([0:Nt jumpTimes1 jumpTimes2]); %IX is the permutation vector
%L = [zeros(1,Nt) L L3(2:end)]; L= L(IX); %sort L over the 'combined' grid
% length of time grid including jump times and discrete times 0,1,...,Nt
lenT = length(T);
% process LT1 over this time grid
LT1 = [zeros(1,Nt+1) L1 zeros(1,length(jumpTimes2)); 0:Nt jumpTimes1 jumpTimes2 ];
LT1 = sortrows(LT1',2)';
% process LT2 over this time grid
LT2 = [zeros(1,Nt+1) L2 zeros(1,length(jumpTimes1)); 0:Nt jumpTimes2 jumpTimes1 ];
LT2 = sortrows(LT2',2)';
% allocate memory for process
E1 = zeros(1,lenT); Y0 = zeros(1,lenT); Y1 = zeros(1,lenT); Y2 = zeros(1,lenT);
% initial values at time 0
E1(1) = mu; Y0(1) = E1(1);
%dt = 1;
for i = 1:lenT-1
    dt = T(i+1) - T(i);
    temp = exp(-1/lambda0*dt);    
    Y0(i+1) = mu*(1-temp) +  temp*Y0(i) + sigma*sqrt((1-temp^2)/2*lambda0)*randn;    
    Y1(i+1) = exp(-1/lambda1*dt)*Y1(i) + LT1(1,i+1);
    Y2(i+1) = exp(-1/lambda2*dt)*Y2(i) + LT2(1,i+1);
    E1(i+1) = Y0(i+1) + negJ(1)*Y1(i+1) + negJ(2)*Y2(i+1);
end

% E1 contains exact simulation of the process, including jump times
% the simulation times are contained in T
% E below contains the observations at times 0, 1, 2, ... the days
%discretise(T,E1);
E = 0; j = 0;
for i = 1:length(T)
    if T(i)==j
        E(j+1) = E1(i);                
        j = j+1;
    end
end

