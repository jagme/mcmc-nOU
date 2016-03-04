% mainMCMC0.m
%
% The following commands can be used to reproduce the results from
%
% 'Bayesian calibration and number of jump components in electricity 
% spot price models' by Jhonny Gonzalez, John Moriarty & Jan Palczewski
%
%
%% ========================================================================
%% Plot deseasonalised data as shown in Figure 1

plotDeasesonalisedData;

%% ========================================================================
%% Fitted seasonal trend coefficients shown in Table A.3

[~,beta1] = fitSeasonalTrend('apxuk1',false);
[~,beta2] = fitSeasonalTrend('eex1',  false);
[~,beta3] = fitSeasonalTrend('apxuk2',false);
[~,beta4] = fitSeasonalTrend('eex2',  false);

fprintf('%10s %12s %12s %12s %10s %13s\n', 'a1','a2','a3','a4','a5','a6')
fprintf('-----------------------------------------------------------------------------------\n')
disp(beta1)
disp(beta2)
disp(beta3)
disp(beta4)

%% ========================================================================
%% Fit 2-OU model to APXUK 2001 - 2006
% Results shown in Table 1, Table 3 and Section 4.2.4
model = '2-OU';
hPar = getHyperparameters(model);  % hyperparameters
x0   = getInitStates(model);       % initial states
n = 2000000;                       % number of mcmc steps
propVar = struct('lambda0',0.05, 'lambda1',0.05); % variance of proposals distributions
NPhi = 5; % number of updates of Phi per mcmc iteration
str = strcat('results/2OU-APXUK1','.mat');  % path to save results
E = fitSeasonalTrend('apxuk1',false);       % deseasonalise the data
% run algorithm
rng(0) % for reproducibility
tic
[par, Lchain, accRate] = gibbs(E, model, n, x0, hPar, propVar, NPhi, str);
toc

% Diagnostics for 2-OU model fit to APXUK 2001 - 2006
% Print posterior mean and SD, and p-values
burnin = 5000; % chain was saved one every 100th iteration
results = diagnostics(str, model, burnin, false);

%% ========================================================================
%% Fit 3-OU model to APXUK 2001 - 2006
% Results shown in Table 2, Table 3 and Section 4.2.4
model = '3-OU';
hPar = getHyperparameters(model);  % hyperparameters
x0   = getInitStates(model);       % initial states
n = 2000000;                       % number of mcmc steps
propVar = struct('lambda0',0.1, 'lambda1',0.05, 'lambda2', 0.05); % variance of proposals distributions
NPhi = 5; % number of updates of Phi per mcmc iteration
str = strcat('results/3OU-APXUK1','.mat');  % path to save results
E = fitSeasonalTrend('apxuk1',false);       % deseasonalise the data
% run algorithm
rng(0) % for reproducibility
tic
[par, accRate] = gibbs3OU(E, model, n, x0, hPar, propVar, NPhi, str);
toc

% Diagnostics for 3-OU model fit to APXUK 2001 - 2006
% Print posterior mean and SD, and p-values
burnin = 5000; % chain was saved one every 100th iteration
results = diagnostics(str, model, burnin, false);

%% ========================================================================
%% Plot Figure 3
% Requires MCMC results for the fitted 2-OU and 3-OU models to APXUK 2001 - 2006
strfile1 = strcat('results/2OU-APXUK1','.mat'); % assumed path to saved results of 2-OU model
strfile2 = strcat('results/3OU-APXUK1','.mat'); % assumed path to saved results of 3-OU model
plotPaths(strfile1, strfile2); % it may take a while to load files if n was large

%% ========================================================================
%% Fit 2-OU model to EEX 2001 - 2006
% Results shown in Table 4, Table 3 and Section 4.3.1
model = '2-OU';
hPar = getHyperparameters(model);  % hyperparameters
x0   = getInitStates(model);       % initial states
n = 2000000;                       % number of mcmc steps
propVar = struct('lambda0',0.05, 'lambda1',0.05); % variance of proposals distributions
NPhi = 5; % number of updates of Phi per mcmc iteration
str = strcat('results/2OU-EEX1','.mat');  % path to save results
E = fitSeasonalTrend('eex1',false);       % deseasonalise the data
% run algorithm
rng(0) % for reproducibility
tic
[par, Lchain, accRate] = gibbs(E, model, n, x0, hPar, propVar, NPhi, str);
toc

% Diagnostics for 2-OU model fit to EEX 2001 - 2006

burnin = 5000; % chain was saved one every 100th iteration
results = diagnostics(str, model, burnin, false);

%% ========================================================================
%% Fit 2-OU-I1 model to EEX 2001 - 2006
% Results shown in Table 4
model = '2-OU-I1';
hPar = getHyperparameters(model);  % hyperparameters
x0   = getInitStates(model);       % initial states
n    = 2000000;                    % number of mcmc steps
propVar = struct('lambda0',0.05, 'lambda1',0.05); % variance of proposals distributions
NPhi = 5; % number of updates of Phi per mcmc iteration
str = strcat('results/2OU-I1-EEX1','.mat');  % path to save results
E = fitSeasonalTrend('eex1',false);       % deseasonalise the data
% run algorithm
rng(0) % for reproducibility
tic
[par, Lchain, accRate] = gibbs(E, model, n, x0, hPar, propVar, NPhi, str);
toc

% Diagnostics for 2-OU-I1 model fit to EEX 2001 - 2006
burnin = 5000; % chain was saved one every 100th iteration
results = diagnostics(str, model, burnin, false);

%% ========================================================================
%% Fit 3-OU-minus model to EEX 2001 - 2006
% Results shown in Table 4,  Table 3 and Section 4.3.1 for 3-OU model
model = '3-OU-';
hPar = getHyperparameters(model);  % hyperparameters
x0   = getInitStates(model);       % initial states
n    = 2000000;                    % number of mcmc steps
propVar = struct('lambda0',0.1, 'lambda1',0.05, 'lambda2', 0.05); % variance of proposals distributions
NPhi = 5; % number of updates of Phi per mcmc iteration
str = strcat('results/3OU-minus-EEX1','.mat');  % path to save results
E = fitSeasonalTrend('eex1',false);       % deseasonalise the data
% run algorithm
rng(0) % for reproducibility
tic
[par, accRate] = gibbs3OU(E, model, n, x0, hPar, propVar, NPhi, str);
toc

% Diagnostics for 3-OU-minus model fit to EEX 2001 - 2006
burnin = 5000; % chain was saved one every 100th iteration
results = diagnostics(str, model, burnin, false);

%% ========================================================================
%% Fit 3-OU-I1 model to EEX 2001 - 2006
% Results shown in Table 4, Table 3, and Section 4.3.1
model = '3-OU-I1';
hPar = getHyperparameters(model);  % hyperparameters
x0   = getInitStates(model);       % initial states
n    = 2000000;                    % number of mcmc steps
propVar = struct('lambda0',0.05, 'lambda1',0.1, 'lambda2', 0.1); % variance of proposals distributions
NPhi = 5; % number of updates of Phi per mcmc iteration
str = strcat('results/3OU-I1-EEX1','.mat');  % path to save results
E = fitSeasonalTrend('eex1',false);       % deseasonalise the data
% run algorithm
rng(0) % for reproducibility
tic
[par, accRate] = gibbs3OU(E, model, n, x0, hPar, propVar, NPhi, str);
toc

% Diagnostics for 3-OU-I1 model fit to EEX 2001 - 2006
burnin = 5000; % chain was saved one every 100th iteration
results = diagnostics(str, model, burnin, false);

%% ========================================================================
%% Plot Figure 2
% Requires MCMC results for the fitted 3-OU model to EEX 2000 - 2006
strfile = strcat('results/3OU-I1-EEX1','.mat');  % assumed file path to results, change as needed
burnin = 5000;
plotJumpSeasonality(strfile, burnin)

%% ========================================================================
%% Prior sensitivity analysis for the 3-OU-I1 model with EEX 2000 - 2006 data
% Results are shown in Table A.4, Prior 2 and Prior 3 columns

%--------------- Estimate model with Prior 2, new prior for sigma 
model = '3-OU-I1';
hPar = getHyperparameters(model);  % hyperparameters
x0   = getInitStates(model);       % initial states
n    = 2000000;                    % number of mcmc steps
propVar = struct('lambda0',0.05, 'lambda1',0.1, 'lambda2', 0.1); % variance of proposals distributions
NPhi = 5; % number of updates of Phi per mcmc iteration
str = strcat('results/3OU-I1-EEX1-sensitivity-sigma','.mat');  % path to save results
E = fitSeasonalTrend('eex1',false);       % deseasonalise the data
% run algorithm
rng(0)
[par, accRate] = gibbs3OU(E, '3-OU-I1-sensitivity-sigma', n, x0, hPar, propVar, NPhi, str);
burnin = 5000;
% display results
results = diagnostics(str, model, burnin, true);

%--------------- % Estimate model with Prior 3, new prior for etai
str = strcat('results/3OU-I1-EEX1-sensitivity-eta','.mat');  % path to save results
rng(0)
[par, accRate] = gibbs3OU(E, '3-OU-I1-sensitivity-eta', n, x0, hPar, propVar, NPhi, str);
% display results
results = diagnostics(str, model, burnin, false);

%% ========================================================================
%% Fit models to 2011 - 2015 data

%% ========================================================================
%% Fit 2-OU model to APXUK 2011 - 2015
% Results shown in Table 5 and 6
model = '2-OU';
hPar = getHyperparameters(model);  % hyperparameters
x0   = getInitStates(model);       % initial states
n    = 2000000;                    % number of mcmc steps
propVar = struct('lambda0',0.05, 'lambda1',0.05); % variance of proposals distributions
NPhi = 5; % number of updates of Phi per mcmc iteration
str = strcat('results/2OU-APXUK2','.mat');  % path to save results
E = fitSeasonalTrend('apxuk2',false);       % deseasonalise the data
% run algorithm
rng(0) % for reproducibility
tic
[par, Lchain, accRate] = gibbs(E, model, n, x0, hPar, propVar, NPhi, str);
toc

% Diagnostics for 2-OU model fit to APXUK 2011 - 2015

burnin = 5000; % chain was saved one every 100th iteration
results = diagnostics(str, model, burnin, false);

%% ========================================================================
%% Fit 2-OU^{-} model to EEX 2011 - 2015
% Results shown in Table 5 and 6
model = '2-OU-';
hPar = getHyperparameters(model);  % hyperparameters
x0   = getInitStates(model);       % initial states
n    = 2000000;                    % number of mcmc steps
propVar = struct('lambda0',0.05, 'lambda1',0.05); % variance of proposals distributions
NPhi = 5; % number of updates of Phi per mcmc iteration
str = strcat('results/2OU-minus-EEX2','.mat');  % path to save results
E = fitSeasonalTrend('eex2',false);       % deseasonalise the data
% run algorithm
rng(0) % for reproducibility
tic
[par, Lchain, accRate] = gibbs(E, model, n, x0, hPar, propVar, NPhi, str);
toc

% Diagnostics for 2-OU^{-} model fit to EEX 2011 - 2015
burnin = 5000; % chain was saved one every 100th iteration
results = diagnostics(str, model, burnin, false);

%% ========================================================================
%% Fit 2-OU-I1 model to EEX 2011 - 2015
% Results shown in Table 5
model = '2-OU-I1';
hPar = getHyperparameters(model);  % hyperparameters
x0   = getInitStates(model);       % initial states
n    = 2000000;                    % number of mcmc steps
propVar = struct('lambda0',0.05, 'lambda1',0.05); % variance of proposals distributions
NPhi = 5; % number of updates of Phi per mcmc iteration
str = strcat('results/2OU-I1-EEX2','.mat');  % path to save results
E = fitSeasonalTrend('eex2',false);       % deseasonalise the data
% run algorithm
rng(0) % for reproducibility
tic
[par, Lchain, accRate] = gibbs(E, model, n, x0, hPar, propVar, NPhi, str);
toc

% Diagnostics for 2-OU-I1 model fit to EEX 2011 - 2015
burnin = 5000; % chain was saved one every 100th iteration
results = diagnostics(str, model, burnin, false);

%% ========================================================================
%% Repeated analysis of the 2-OU model on simulated data
% Results shown in Table A.2
strfile = 'results/2c-estimatesSimT1000.mat'; % file path to save results
burnin = 3000;

%{
% the true simulation values are
mu      = 1;    % level of mean reversion
lambda0 = 8;    % time to mean reversion
sigma   = 0.1;  % volatility
lambda1 = 2;    % time to mean reversion 
beta    = 0.7;  % mean jump size

eta varies over [0.05 0.1 0.2 0.3 ]
%}
% this function may take long to complete, > 15 hours
% use the parallel toolbox if available to improve computational time
repeatAnalysis(strfile);
results = processRepeatAnalysis(strfile, burnin);

%% ========================================================================
%% Compare ACF with different number of updates of Phi per MCMC iteration
% as shown in Figure A.1 
% -------------- simulate 3-OU process
w = [1 1];  % sign of each jump component
Nt = 1000;
mu      = 1;    sigma   = 0.15;
lambda0 = 8;    lambda1 = 3; lambda2 = 0.5;
eta1    = 0.1;  eta2    = 0.05;
beta1   = 0.5;  beta2   = 1;
rng(0)
E = simulate3OUModel(Nt,mu, lambda0, sigma, lambda1,lambda2, eta1, eta2, beta1, beta2, w);

model = '3-OU';
hPar = getHyperparameters(model);  % hyperparameters
x0   = getInitStates(model);       % initial states
n = 1000000;                       % number of mcmc steps
propVar = struct('lambda0',0.1, 'lambda1',0.05, 'lambda2', 0.05); % variance of proposals distributions

% -------------- estimate 3-OU model with 1 update of Phi per MCMC iter.
rng(0)
strfile1 = strcat('results/3OU-simulationNPhi_1','.mat');  % path to save results
NPhi = 1; % number of updates of Phi per mcmc iteration
[par1, accRate] = gibbs3OU(E, model, n, x0, hPar, propVar, NPhi, strfile1);

% -------------- estimate 3-OU model with 5 updates of Phi per MCMC iter.
rng(0)
strfile2 = strcat('results/3OU-simulationNPhi_5','.mat');  % path to save results
NPhi = 5; % number of updates of Phi per mcmc iteration
[par2, accRate] = gibbs3OU(E, model, n, x0, hPar, propVar, NPhi, strfile2);

% -------------- plot ACFs of eta
burnin = 5000;
figure;
hold on
ACF1 = autocorr(par1(burnin:end,5),800);
ACF2 = autocorr(par2(burnin:end,5),800);
plot(ACF1)
plot(ACF2)
grid on, box off, ylabel('ACF of \eta_1')
legend('One update of \Phi per MCMC iteration','Five updates of \Phi per MCMC iteration')
h = gca; h.FontSize = 11;

%% ========================================================================
%% Summary of process Phi = (Phi1, Phi2) of 3-OU model
% as shown in Figure A.2

% -------------- simulate 3-OU process and save true processes Phi1, Phi2
w = [1 1];      % sign of each jump component
Nt = 1000;
mu      = 1;    sigma   = 0.15;
lambda0 = 8;    lambda1 = 3;    lambda2 = 0.5;
eta1    = 0.1;  eta2    = .05;
beta1   = 0.5;  beta2   = 1;
rng(0)
[E,LT1,LT2] = simulate3OUModel(Nt,mu, lambda0, sigma, lambda1,lambda2, eta1, eta2, beta1, beta2, w);
strfile1 = 'results/3OU-simulationLT1LT2.mat';
save(strfile1,'E','LT1','LT2')

% -------------- estimate 3-OU model using the corresponding observed process
model = '3-OU';
hPar = getHyperparameters(model);  % hyperparameters
x0   = getInitStates(model);       % initial states
n    = 1000000;                    % number of mcmc steps
NPhi = 5; % number of updates of Phi per mcmc iteration
propVar = struct('lambda0',0.1, 'lambda1',0.05, 'lambda2', 0.05); % variance of proposals distributions
strfile2 = strcat('results/3OU-simulationrng0','.mat');  % path to save results
% run algorithm
rng(0)
tic
[par, accRate] = gibbs3OU(E, model, n, x0, hPar, propVar, NPhi, strfile2);
toc

% -------------- plot true Phi vs estimated Phi
plotPhi(strfile1,strfile2)

