function parM = changePar(Nt, eta, Niterations)
%CHANGEPAR Fit the 2-OU model to simulated data using pre-specified 
% parameters mu, lambda0, sigma, lambda1, beta1 &
% eta1 passed as input argument. 
%
% parM = changePar(Nt, eta, Niterations)
% 
% Input arguments:
%
%   Nt - Number of timesteps for simulated path
%
%   eta - parameter eta: jump intensity rate of L
%
%   Niterations - number of MCMC draws
%
% Output arguments:
%
%   parM - chain of parameters
%

model = '2-OU';
w1 = 1;
% set simulation values
mu = 1;         % level of mean reversion
lambda0 = 8;    % time to mean reversion
sigma = 0.1;    % volatility
lambda1 = 2;    % time to mean reversion 
jSizePar = 0.7; % mean jump size

E = simulateModel(Nt, mu, lambda0, sigma, lambda1, eta, jSizePar, w1);

hPar = getHyperparameters(model);  % hyperparameters
hPar.Beta = 1/eta;
x0   = getInitStates(model);       % initial states

propVar = struct('lambda0',0.05, 'lambda1',0.05); % variance of proposals distributions
NPhi = 1; % number of updates of Phi per mcmc iteration
str = ''; % do not save MCMC results from gibbs.m

parM = gibbs(E, model, Niterations, x0, hPar, propVar, NPhi, str);
