function L = likelihoodPP(Tmax,Npp,tau,t0,period,d,theta)
%LIKELIHOODPP Log-likelihood of inhomogeneous Poisson process with jump
%intensity given in file intensityFun.m
%
% L = likelihoodPP(Tmax,Npp,tau,t0,period,d,theta)
% 
% Input arguments:
%
%   Tmax    - Process lives on [0, Tmax]
%   Npp     - Number of jumps in [0, Tmax]
%   tau     - jump times of the realisation of the Poisson process

%  parameters to intensity function:
%
%   t0      - phase
%   period  - period of intensity function
%   d       - controls dispersion of jumps around peaking levels
%   theta   - max. expected number of jumps

% handle for intensity rate function
fhandle = @(t)intensityFun(t, t0, period, d, theta);
% compute integral on interval [0,Tmax]
rate = integral(fhandle,0,Tmax);

% to avoid overflow, don't compute likelihood here but only return
% expression inside the exponential of the likelihood ratio

L = sum(log(intensityFun(tau, t0, period, d, theta))) - rate;