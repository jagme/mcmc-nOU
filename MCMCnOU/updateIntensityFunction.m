function [delta1, eta1, theta1, accRate] = updateIntensityFunction(Tmax,Tau,t0,period,delta,eta)
%UPDATEINTENSITYFUNCTION MH step for parameters of the jump intensity
%function in intensityFun.m
%
% [delta1, eta1, theta1, accRate] = updateIntensityFunction(Tmax,Tau,t0,period,delta,eta)
%
% Input arguments:
%
%   Tmax - Poisson process lives on [0, Tmax]
%
%   Tau - observations of the jump times of the Poisson process
%
%   period - period of intensity function
%
%   t0, delta, eta -  current states of parameters
%
% Output arguments:
%
%   delta1, eta1, theta1 - new states of the parameters delta, eta and theta
%
%   accRate - logical vector indicating whether the proposals for each
%   parameter where accepted
%
N = 2;
% initial states
chain = NaN(N,3);
chain(1,:) = [delta eta t0]; % delta1, eta1, theta1
% hyperparameters
Ad = 1;
sigmaRW = 0.2;
Npp = length(Tau);
accRate = zeros(3,1); % acceptance rates
% loop MH
for i = 2:N
    % delta
    x = chain(i-1,1);
    y = exp(normrnd(log(x),sigmaRW));    
    ratex = likelihoodPP(Tmax,Npp,Tau,t0,period,x,chain(i-1,2));
    ratey = likelihoodPP(Tmax,Npp,Tau,t0,period,y,chain(i-1,2));
    L =  exp(ratey-ratex) *   exppdf(y,Ad)/exppdf(x,Ad) * y/x;
    if rand <= L
        chain(i,1) = y;
        accRate(1) = accRate(1) + 1;
    else
        chain(i,1) = x;
    end     
    
    % eta
    x = chain(i-1,2);
    y = exp(normrnd(log(x),sigmaRW));    
    ratex = likelihoodPP(Tmax,Npp,Tau,t0,period,chain(i-1,1),x);
    ratey = likelihoodPP(Tmax,Npp,Tau,t0,period,chain(i-1,1),y);    
    L =  exp(ratey-ratex) *    exppdf(y,Ad)/exppdf(x,Ad) * y/x;
    if rand <= L
        chain(i,2) = y; 
        accRate(2) = accRate(2) + 1;
    else
        chain(i,2) = x;
    end
    
    % theta
    x = chain(i-1,3);
    y = exp(normrnd(log(x),sigmaRW));
    
    ratex = likelihoodPP(Tmax,Npp,Tau,x,period,chain(i-1,1),chain(i-1,2));
    ratey = likelihoodPP(Tmax,Npp,Tau,y,period,chain(i-1,1),chain(i-1,2));
    L =  exp(ratey-ratex) *   y/x *and(y>65,y<195)*and(x>65,x<195); 
    % note the peaking level multiple (260/2, or 6-month here) should live
    % either on a suitable interval e.g. [0,130] or [130, 260], or [65, 195] etc. 
    % But the pair x = 131, y = 132 has greater likelihood than x = 1, y = 2, 
    % although they are the same angles (up to period = 130 days). 
    % This is because of the ratio y/x. Therefore we use the interval
    % [65,195]. There is no identifiability issue here.
    if rand <= L
        chain(i,3) = y; 
        accRate(3) = accRate(3) + 1;
    else
        chain(i,3) = x;
    end
   
end

delta1 = chain(N,1);    % delta1, eta1, theta1
eta1   = chain(N,2);
theta1 = chain(N,3);
    
