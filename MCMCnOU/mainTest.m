% test file for MCMC algorithms with simulated data

% The goal is to check that the true simulation values are recovered by 
% the algorithms. To get more posterior samples and/or study ACFs
% increase the number of mcmc steps n to n > 1000000

% models
% 2-OU, 2-OU-, 2-OU-I1, 3-OU, 3-OU-, 3-OUI1

% 2-OU      one positive jump component
% 2-OU-     one negative jump component
% 2-OU-I1   one positive jump component, L1 is time inhomogeneous

% 3-OU      two positive jump components
% 3-OU-     Y1 positive, Y2 negative
% 3-OU-I1   Y1 positive, Y2 negative, L1 is time inhomogeneous

%% ========================================================================
%% ===================== test 2-OU models =================================
%% test 2-OU model
rng(1)
model = '2-OU';
w  = 1;         % sign of jump component Y1
Nt = 1000;      % number of timesteps
% set simulation values
mu      = 1;    % level of mean reversion
lambda0 = 8;    % time to mean reversion
sigma   = 0.1;  % volatility
lambda1 = 2;    % time to mean reversion 
eta     = 0.05; % jump rate 
beta    = 0.7;  % mean jump size

[E, T, L] = simulateModel(Nt, mu, lambda0, sigma, lambda1, eta, beta, w);

n = 100000;                        % number of mcmc steps
hPar = getHyperparameters(model);  % hyperparameters
x0 = getInitStates(model);         % initial states
propVar = struct('lambda0',0.05, 'lambda1',0.05); % variance of proposals distributions
str = strcat('results/2OU-sim','.mat'); % path to save results
NPhi = 1; % number of updates of Phi per mcmc iteration
% run algorithm
tic
[par, Lchain, accRate] = gibbs(E, model, n, x0, hPar, propVar, NPhi, str);
toc
burnin = floor(n/100*0.2); % chain was saved one every 100th iteration
results = diagnostics(str, model, burnin, true);

%% test 2-OU- model
rng(1)
model = '2-OU-';
w  = -1;        % sign of jump component Y1
Nt = 1000;      % number of timesteps
% set simulation values
mu      = 1;    % level of mean reversion
lambda0 = 8;    % time to mean reversion
sigma   = 0.1;  % volatility
lambda1 = 2;    % time to mean reversion 
eta     = 0.05; % jump rate 
beta    = 0.7;  % mean jump size

[E, T, L] = simulateModel(Nt, mu, lambda0, sigma, lambda1, eta, beta, w);

n = 100000;                        % number of mcmc steps
hPar = getHyperparameters(model);  % hyperparameters
x0 = getInitStates(model);         % initial states
propVar = struct('lambda0',0.05, 'lambda1',0.05); % variance of proposals distributions
str = strcat('results/2OUminus-sim','.mat'); % path to save results
NPhi = 1; % number of updates of Phi per mcmc iteration
% run algorithm
[par, Lchain, accRate] = gibbs(E, model, n, x0, hPar, propVar, NPhi, str);

burnin = floor(n/100*0.2); % chain was saved one every 100th iteration
results = diagnostics(str, model, burnin, true);

%% test 2-OU-I1 model
rng(1)
model = '2-OU-I1';
w  = 1;             % sign of jump component Y1
Nt = 1000;          % number of timesteps
% set simulation values
mu      = 1;        % level of mean reversion
lambda0 = 8;        % time to mean reversion
sigma   = 0.1;      % volatility
lambda1 = 2;        % time to mean reversion 
beta    = 0.7;      % mean jump size
% parameters for jump intensity function I1
I1.t0     = 120;    % theta
I1.delta  = 0.2;    % delta
I1.eta    = 0.5;    % eta
I1.period = 260;    % period 

[E, T, L] = simulateModel(Nt, mu, lambda0, sigma, lambda1, I1, beta, w);

n = 200000;                        % number of mcmc steps
hPar = getHyperparameters(model);  % hyperparameters
x0 = getInitStates(model);         % initial states
propVar = struct('lambda0',0.05, 'lambda1',0.05); % variance of proposals distributions
str = strcat('results/2OUI1-sim','.mat'); % path to save results
NPhi = 1; % number of updates of Phi per mcmc iteration
% run algorithm
tic
[par, Lchain, accRate] = gibbs(E, model, n, x0, hPar, propVar, NPhi, str);
toc
burnin = floor(n/100*0.2); % chain was saved one every 100th iteration
results = diagnostics(str, model, burnin, true);

%% ========================================================================
%% ===================== test 3-OU models =================================
%% test 3-OU model
% -------------- simulate 3-OU process
model = '3-OU';
w     = [1 1]; % sign of each jump component
Nt = 1000;
mu      = 1;    sigma   = 0.15;
lambda0 = 8;    lambda1 = 3; lambda2 = 0.5;
eta1    = 0.1;  eta2    = 0.05;
beta1   = 0.5;  beta2   = 1;

rng(0)
E = simulate3OUModel(Nt,mu, lambda0, sigma, lambda1,lambda2, eta1, eta2, beta1, beta2, w);

hPar = getHyperparameters(model);  % hyperparameters
x0   = getInitStates(model);       % initial states
n = 100000;                        % number of mcmc steps
propVar = struct('lambda0',0.1, 'lambda1',0.05, 'lambda2', 0.05); % variance of proposals distributions

% -------------- estimate 3-OU model with 1 update of Phi per MCMC iter.
rng(0)
str = strcat('results/3OU-simulationNPhi_1','.mat');  % path to save results
NPhi = 1; % number of updates of Phi per mcmc iteration
tic
[par1, accRate] = gibbs3OU(E, model, n, x0, hPar, propVar, NPhi, str);
toc
burnin = floor(n/100*0.3); % chain was saved one every 100th iteration
results = diagnostics(str, model, burnin, true);

%% test 3-OU- model
% -------------- simulate 3-OU- process
model = '3-OU-';
w     = [1 -1]; % sign of each jump component
Nt = 1000;
mu      = 1;    sigma   = 0.15;
lambda0 = 8;    lambda1 = 3; lambda2 = 0.5;
eta1    = 0.1;  eta2    = 0.05;
beta1   = 0.5;  beta2   = 1;

rng(0)
E = simulate3OUModel(Nt,mu, lambda0, sigma, lambda1,lambda2, eta1, eta2, beta1, beta2, w);

hPar = getHyperparameters(model);  % hyperparameters
x0   = getInitStates(model);       % initial states
n = 100000;                        % number of mcmc steps
propVar = struct('lambda0',0.1, 'lambda1',0.05, 'lambda2', 0.05); % variance of proposals distributions

% -------------- estimate 3-OU model with 1 update of Phi per MCMC iter.
rng(0)
str = strcat('results/3OUminus-simulationNPhi_1','.mat');  % path to save results
NPhi = 1; % number of updates of Phi per mcmc iteration
[par1, accRate] = gibbs3OU(E, model, n, x0, hPar, propVar, NPhi, str);

burnin = floor(n/100*0.3); % chain was saved one every 100th iteration
results = diagnostics(str, model, burnin, true);

%% test 3-OU-I1 model
% -------------- simulate 3-OU-I1 process
model = '3-OU-I1';
w     = [1 -1]; % sign of each jump component
Nt = 1000;
mu      = 1;    sigma   = 0.15;
lambda0 = 8;    lambda1 = 3; lambda2 = 0.5;
% eta1    = 0.1;  
I1.t0     = 120; % parameters of I1, jump intensity function
I1.eta    = 0.3;
I1.delta  = 0.1;
I1.period = 260;
eta2      = 0.05;
beta1     = 0.5;  beta2   = 1;

rng(0)
E = simulate3OUModel(Nt,mu, lambda0, sigma, lambda1,lambda2, I1, eta2, beta1, beta2, w);

hPar = getHyperparameters(model);  % hyperparameters
x0   = getInitStates(model);       % initial states
n = 100000;                        % number of mcmc steps
propVar = struct('lambda0',0.1, 'lambda1',0.05, 'lambda2', 0.05); % variance of proposals distributions

% -------------- estimate 3-OU model with 1 update of Phi per MCMC iter.
rng(0)
str = strcat('results/3OUI1-simulationNPhi_1','.mat');  % path to save results
NPhi = 1; % number of updates of Phi per mcmc iteration
[par1, accRate] = gibbs3OU(E, model, n, x0, hPar, propVar, NPhi, str);

burnin = floor(n/100*0.3); % chain was saved one every 100th iteration
results = diagnostics(str, model, burnin, true);
