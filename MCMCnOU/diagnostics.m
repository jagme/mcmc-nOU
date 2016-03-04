function results = diagnostics(file, model, burnin, plotTr)
%DIAGNOSTICS Compute posterior mean, posterior standard deviation, and
%Bayesian p-values, and optionally plot trace, posterior and ACFs for each
%parameter
%
% Input arguments:

% file - file path to saved mcmc results from either gibbs.m or gibbs3OU.m
%
% model - string specifying the fitted model
%       '2-OU', '2-OU-', '2-OU-I1'
%       '3-OU', '3-OU-', '3-OU-I1'
%
% burnin - the burn-in period, the first burnin samples will be removed to
% compute the diagnostics
%
% plotTr - boolean variable, if true traces, posteriors and ACFs are
% plotted

switch model
    case {'2-OU', '2-OU-', '2-OU-I1'}
        results = d2OU(file, model, burnin, plotTr);                
        
    case {'3-OU', '3-OU-','3-OU-I1'}
        results = d3OU(file, model, burnin, plotTr);                
        
    otherwise
        error('Model has been mis-specified')
end

end % end diagnostics

% diagnostics for 2-ou models
function results = d2OU(file, model, burnin, plotTr)

load(file)
Tmax = max(T); % max observation time

% posterior mean and SD, optinally plot trace, posterior density and ACF
switch model
    case {'2-OU', '2-OU-'}
        post_mean = stats2OU(par, E, burnin, plotTr);        
    case '2-OU-I1'
        post_mean = stats2OUI1(par, E, burnin, plotTr);
end

% compute p-values
pvWiener = dWiener2OU(E, par, Lchain, burnin, model);
fprintf('P-value for Y0: %f\n',pvWiener)
pvJump = dJump2OU(par, Lchain, Tmax, burnin, model);
fprintf('P-value for jump times of Y1: %f\n',pvJump(1))
fprintf('P-value for jump sizes of Y1: %f\n',pvJump(2))
results.postStats = post_mean;
results.pvWiener  = pvWiener;
results.pvJump    = pvJump;

end % end d2OU

% diagnostics for 3-ou models
function results = d3OU(file, model, burnin, plotTr)

load(file)
Tmax = max(T); % max observation time

% posterior mean and SD, optinally plot trace, posterior density and ACF
switch model
    case {'3-OU', '3-OU-'}
        post_mean = stats3OU(par, E, burnin, plotTr);
    case '3-OU-I1'
        post_mean = stats3OUI1(par, E, burnin, plotTr);
    
end
% compute p-values
pvWiener = dWiener3OU(E, par, Lchain_1, Lchain_2, burnin, model);
fprintf('P-value for Y0: %f\n',pvWiener)
pvJump = dJump3OU(par, Lchain_1, Lchain_2, Tmax, burnin, model);
fprintf('P-value for jump times of Y1: %f\n',pvJump(1))
fprintf('P-value for jump sizes of Y1: %f\n',pvJump(2))
fprintf('P-value for jump times of Y2: %f\n',pvJump(3))
fprintf('P-value for jump sizes of Y2: %f\n',pvJump(4))

results.postStats = post_mean;
results.pvWiener  = pvWiener;
results.pvJump    = pvJump;

end % end d3OU

% compute statistics (mainly posterior mean and std) of chain
% optionally plot trace, posterior and ACF
function postProperties = stats2OU(par, E, burnin, plotTr)

if burnin >= length(par)
   error('Burn-in period is too long, must be smaller than %d',length(par)) 
end

par = par(burnin:end,:); % remove burn in period

Tmax     = length(E) - 1;

mu      = par(:,1);
lambda0 = par(:,2);
sigma   = par(:,3);
lambda1 = par(:,4);
eta     = par(:,5)/Tmax;
beta    = par(:,6);

%         [1    2         3      4         5      6         7     8   ]
strpar = {'mu' 'sigma^2' 'rho0' 'lambda0' 'rho1' 'lambda1' 'eta' 'beta'};
par2   = [mu sigma.^2 exp(-1./lambda0) lambda0 exp(-1./lambda1) lambda1 eta beta];
postProperties = [mean(par2)' std(par2)'];

format shortg
fprintf('Posterior properties\n')
fprintf('%10s \t %s \t %s \n', 'Parameter', 'Post. mean', 'Post. SD')
fprintf('------------------------------------------\n')
for i = 1:length(postProperties)
    fprintf('%10s \t %f \t %f \n', strpar{i}, postProperties(i,1), postProperties(i,2))
end
fprintf('------------------------------------------\n')

% plot trace, posterior and ACF?
if plotTr
    figure; H = 3; V = 6; burn = 1;
    plotR(mu,       burn, [V H 1],   '\mu')
    plotR(lambda0,  burn, [V H 4],   '\lambda_0')
    plotR(sigma.^2, burn, [V H 7],   '\sigma^2')
    plotR(lambda1,  burn, [V H 10],  '\lambda_1')
    plotR(eta,      burn, [V H 13],  '\eta')
    plotR(beta,     burn, [V H 16],  '\beta')
    
    drawnow
end

end % end stats2OU

function postProperties = stats2OUI1(par, E, burnin, plotTr)

if burnin >= length(par)
   error('Burn-in period is too long, must be smaller than %d',length(par)) 
end

par = par(burnin:end,:); % remove burn in period

Tmax     = length(E) - 1;

mu      = par(:,1);
lambda0 = par(:,2);
sigma   = par(:,3);
lambda1 = par(:,4);
beta    = par(:,5);
delta1  = par(:,6);
eta     = par(:,7);
theta1  = par(:,8);


%         [1    2         3      4         5      6         7     8   ]
strpar = {'mu' 'sigma^2' 'rho0' 'lambda0' 'rho1' 'lambda1' 'eta' 'beta' 'theta' 'delta'};
par2   = [mu sigma.^2 exp(-1./lambda0) lambda0 exp(-1./lambda1) lambda1 eta beta theta1 delta1];
postProperties = [mean(par2)' std(par2)'];

format shortg
fprintf('Posterior properties\n')
fprintf('%10s \t %s \t %s \n', 'Parameter', 'Post. mean', 'Post. SD')
fprintf('------------------------------------------\n')
for i = 1:length(postProperties)
    fprintf('%10s \t %f \t %f \n', strpar{i}, postProperties(i,1), postProperties(i,2))
end
fprintf('------------------------------------------\n')

% plot trace, posterior and ACF?
if plotTr
    figure; H = 3; V = 4; burn = 1;
    plotR(mu,       burn, [V H 1],   '\mu')
    plotR(lambda0,  burn, [V H 4],   '\lambda_0')
    plotR(sigma.^2, burn, [V H 7],   '\sigma^2')
    plotR(lambda1,  burn, [V H 10],  '\lambda_1')
    
    figure; H = 3; V = 4; burn = 1;
    plotR(eta,      burn, [V H 1],  '\eta_1')
    plotR(beta,     burn, [V H 4],  '\beta_1')
    plotR(theta1,   burn, [V H 7],  '\theta_1')
    plotR(delta1,   burn, [V H 10], '\delta_1')
    
    drawnow
end

end % end stats2OUI1


% compute statistics (mainly posterior mean and std) of chain
% optional plot trace, posterior and ACF
function postProperties = stats3OU(par, E, burnin, plotTr)

if burnin >= length(par)
   error('Burn-in period is too long, must be smaller than %d',length(par)) 
end

par = par(burnin:end,:); % remove burn in period

Tmax     = length(E) - 1;

mu      = par(:,1);
lambda0 = par(:,2);
sigma   = par(:,3);
lambda1 = par(:,4);
eta1     = par(:,5)/Tmax;
lambda2 = par(:,6);
eta2    = par(:,7)/Tmax;
beta1   = par(:,8);
beta2   = par(:,9);

strpar = {'mu' 'sigma^2' 'rho0' 'lambda0' 'rho1' 'lambda1' 'rho2' 'lambda2' 'eta1' 'eta2' 'beta1' 'beta2'};

par2   = [mu sigma.^2 exp(-1./lambda0) lambda0 exp(-1./lambda1) lambda1...
    exp(-1./lambda2) lambda2 eta1 eta2 beta1 beta2];
postProperties = [mean(par2)' std(par2)'];

format shortg
fprintf('Posterior properties\n')
fprintf('%10s \t %s \t %s \n', 'Parameter', 'Post. mean', 'Post. SD')
fprintf('------------------------------------------\n')
for i = 1:length(postProperties)
    fprintf('%10s \t %f \t %f \n', strpar{i}, postProperties(i,1), postProperties(i,2))
end
fprintf('------------------------------------------\n')
% plot trace, posterior and ACF?
if plotTr
    figure; H = 3; V = 4; burn = 1;
    plotR(mu,       burn, [V H 1],   '\mu')
    plotR(lambda0,  burn, [V H 4],   '\lambda_0')
    plotR(sigma.^2, burn, [V H 7],   '\sigma^2')
    plotR(lambda1,  burn, [V H 10],  '\lambda_1')
    
    figure; H = 3; V = 5;
    plotR(lambda2,   burn, [V H 1],  '\lambda_2')
    plotR(eta1,      burn, [V H 4],  '\eta_1')
    plotR(eta2,      burn, [V H 7],  '\eta_2')
    plotR(beta1,     burn, [V H 10],  '\beta_1')
    plotR(beta1,     burn, [V H 13],  '\beta_2')
    
    drawnow
end

end % end stats3OU

% compute statistics (mainly posterior mean and std) of chain
% optional plot trace, posterior and ACF
function postProperties = stats3OUI1(par, E, burnin, plotTr)

if burnin >= length(par)
   error('Burn-in period is too long, must be smaller than %d',length(par)) 
end

par = par(burnin:end,:); % remove burn in period

Tmax     = length(E) - 1;

mu      = par(:,1);
lambda0 = par(:,2);
sigma   = par(:,3);
lambda1 = par(:,4);
eta1     = par(:,5);
lambda2 = par(:,6);
eta2    = par(:,7)/Tmax;
beta1   = par(:,8);
beta2   = par(:,9);
delta1  = par(:,10);
t0      = par(:,11);

strpar = {'mu' 'sigma^2' 'rho0' 'lambda0' 'rho1' 'lambda1' 'rho2' 'lambda2' 'eta1' 'eta2' 'beta1' 'beta2','theta1','delta1'};

par2   = [mu sigma.^2 exp(-1./lambda0) lambda0 exp(-1./lambda1) lambda1...
    exp(-1./lambda2) lambda2 eta1 eta2 beta1 beta2 t0 delta1];
postProperties = [mean(par2)' std(par2)'];

format shortg
fprintf('Posterior properties\n')
fprintf('%10s \t %s \t %s \n', 'Parameter', 'Post. mean', 'Post. SD')
fprintf('------------------------------------------\n')
for i = 1:length(postProperties)
    fprintf('%10s \t %f \t %f \n', strpar{i}, postProperties(i,1), postProperties(i,2))
end
fprintf('------------------------------------------\n')

% plot trace, posterior and ACF?
if plotTr
    figure; H = 3; V = 4; burn = 1;
    plotR(mu,       burn, [V H 1],   '\mu')
    plotR(lambda0,  burn, [V H 4],   '\lambda_0')
    plotR(sigma.^2, burn, [V H 7],   '\sigma^2')
    plotR(lambda1,  burn, [V H 10],  '\lambda_1')
    
    figure; H = 3; V = 4;
    plotR(lambda2,   burn, [V H 1],  '\lambda_2')
    plotR(eta1,      burn, [V H 4],  '\eta_1')
    plotR(eta2,      burn, [V H 7],  '\eta_2')
    plotR(beta1,     burn, [V H 10],  '\beta_1')
    
    figure; H = 3; V = 3;
    plotR(beta2,     burn, [V H 1],  '\beta_2')
    plotR(t0,        burn, [V H 4],  '\theta_1')
    plotR(delta1,    burn, [V H 7],  '\delta_1')
    
    drawnow
end

end % end stats3OUI1

% compute p-value for residuals of implied Wiener process
function meanpVal = dWiener2OU(E, par, Lchain, burnin,model)

Nt   = length(E);
Nmax = length(Lchain) - burnin;
step = 5;
pVal = NaN(Nmax/step,1);

switch model
    case '2-OU'
        w1 = 1;
    case '2-OU-I1'
        w1 = 1;
    case '2-OU-'
        w1 = -1;
    otherwise
        error('Model has been mis-specified')
end

for t = 1:Nmax/step
    % get induced OU process
    t2 = t*step + burnin;
    Y2 = getY2(par(t2,4),Lchain{t2},Nt);
    Y1 = E - w1*Y2;
    % get residuals of Wiener process
    mu = par(t2,1);
    sigma = par(t2,3);
    lambda0 = par(t2,2);
    res = getWResiduals(Y1, mu, lambda0, sigma);
    
    % One-sample Kolmogorov-Smirnov test
    % h = 1 indicate rejection of the null of no autocorrelation in favor of the alternative.
    % h = 0 indicate a failure to reject the null, evidence there's autocorrelation
    [h, pValue] = kstest(res);
    pVal(t)    = pValue;
    
end

meanMedian = [ mean(pVal) median(pVal)];
meanpVal = meanMedian(1);

end % end dWiener2OU

% compute p-value for residuals of implied Wiener process
function meanpVal = dWiener3OU(E, par, Lchain_1, Lchain_2, burnin, model)

Nt   = length(E);
Nmax = length(Lchain_1) - burnin;
step = 5;
pVal = NaN(Nmax/step,1);

switch model
    case '3-OU'
        w1 = 1; w2 = 1;
    case '3-OU-I1'
        w1 = 1; w2 = -1;
    case '3-OU-'
        w1 = 1; w2 = -1;
    otherwise
        error('Model has been mis-specified')
end

for t = 1:Nmax/step
    % get induced OU process
    t2 = t*step + burnin;
    Y2 = getY2(par(t2,4),Lchain_1{t2},Nt);
    Y3 = getY2(par(t2,6),Lchain_2{t2},Nt);
    Y1 = E - w1*Y2 - w2*Y3;           % change here for positive and negative jumps
    % get residuals of Wiener process
    mu = par(t2,1);
    sigma = par(t2,3);
    lambda0 = par(t2,2);
    res = getWResiduals(Y1, mu, lambda0, sigma);
    
    % One-sample Kolmogorov-Smirnov test
    % h = 1 indicate rejection of the null of no autocorrelation in favor of the alternative.
    % h = 0 indicate a failure to reject the null, evidence there's autocorrelation
    [h, pValue] = kstest(res);
    pVal(t)    = pValue;
    
end

meanMedian = [ mean(pVal) median(pVal)];
meanpVal = meanMedian(1);

end % end dWiener3OU

% obtain the residuals of the implied Gaussian OU process Y
function res = getWResiduals(Y, mu, lambda0, sigma)

Nt = length(Y);
W  = NaN(1, Nt - 1);
temp = exp(-1/lambda0);
temp2 = sqrt(sigma^2*lambda0/2*(1-exp(-2/lambda0)));

for i = 1:Nt-1
    W(i) = (Y(i+1) - mu - (Y(i) - mu)*temp)/temp2;
end

res = W;

end

function pv = dJump2OU(par, Lchain, Tmax, burnin, model)

Nmax = length(Lchain) - burnin;
step = 1;
pVal = NaN(2,Nmax/step,1);

alpha = 0.1; % Significance level for test

switch model
    case {'2-OU', '2-OU-'}
        for t = 1:Nmax/step
            t2 = t*step + burnin;      % index
            testData  = exprnd(par(t2,6),10000,1);        % for jump sizes
            testData2 = exprnd( Tmax/par(t2,5) ,10000,1); % for jump times
            
            Phi = Lchain{t2};      % Phi at iteration t2
            % get jump sizes
            jSizes = Phi(1,:);
            idx = jSizes>0;
            jSizes = jSizes(idx);  % the jump sizes
            % get jump times
            jTimes = Phi(2,idx);
            jTimes = diff(jTimes); % the interarrival times
            
            if numel(jTimes) > 1
                [h,p] = kstest2(testData2,jTimes,'Alpha',alpha); % test jump times
                pVal(1,t) = p;
                [h,p] = kstest2(testData,jSizes,'Alpha',alpha); % test jump sizes
                pVal(2,t) = p;
            else
                pVal(1,t) = NaN;
                pVal(2,t) = NaN;
            end
            
        end
        
    case '2-OU-I1'
        
        for t = 1:Nmax/step
            t2 = t*step + burnin;      % index
            testData = exprnd(par(t2,5),10000,1);   % for jump sizes
            % test data from inhomogeneous CPP process
            [~,testData2] = compoundPoisson(Tmax*5,par(t2,8),260,par(t2,6),par(t2,7),par(t2,5),2);
            testData2 = mod(testData2, 260); % yearly, 260 days
            
            Phi = Lchain{t2};    % Phi iteration
            % get jump sizes
            jSizes = Phi(1,:);
            idx = jSizes>0;
            jSizes = jSizes(idx);
            % get jump times
            jTimes = Phi(2,idx);
            jTimes = mod((jTimes), 260);
            
            if numel(jTimes)>1
                [h,p] = kstest2(testData2,jTimes,'Alpha',alpha); % test jump times
                pVal(1,t) = p;
                [h,p] = kstest2(testData,jSizes,'Alpha',alpha);  % test jump sizes
                pVal(2,t) = p;
            else
                pVal(1,t) = NaN;
                pVal(2,t) = NaN;
            end
            
        end
end

pv = nanmean(pVal, 2); % compute mean p-value


end % end dJump2OU

function pv = dJump3OU(par, Lchain_1, Lchain_2, Tmax, burnin, model)

Nmax = length(Lchain_1) - burnin;
step = 1;
pVal = NaN(4,Nmax/step,1);

alpha = 0.1; % Significance level for test

switch model
    case {'3-OU', '3-OU-'}
        for t = 1:Nmax/step
            t2 = t*step + burnin;     % index
            testData  = exprnd(par(t2,8),10000,1);
            testData2 = exprnd( Tmax/par(t2,5) ,10000,1);            
            testData3 = exprnd(par(t2,9),10000,1);
            testData4 = exprnd( Tmax/par(t2,7) ,10000,1);
            % ----------- Lchain1
            Phi = Lchain_1{t2};     % Phi iteration
            % get jump sizes
            jSizes = Phi(1,:);
            idx = jSizes>0;
            jSizes = jSizes(idx);
            % get jump times
            jTimes = Phi(2,idx);
            jTimes = diff(jTimes);            
            
            if numel(jTimes) > 1                                                
                [h,p] = kstest2(testData2,jTimes,'Alpha',alpha); % test jump times                
                pVal(1,t) = p;
                [h,p] = kstest2(testData,jSizes,'Alpha',alpha);  % test jump sizes
                pVal(2,t) = p;
            else
                pVal(1,t) = NaN;
                pVal(2,t) = NaN;
            end
            % ----------- Lchain2
            Phi = Lchain_2{t2};     % Phi iteration
            % get jump sizes
            jSizes = Phi(1,:);
            idx = jSizes>0;
            jSizes = jSizes(idx);
            % get jump times
            jTimes = Phi(2,idx);
            jTimes = diff(jTimes);
            if numel(jTimes) > 1                
                [h,p] = kstest2(testData4,jTimes,'Alpha',alpha); % test jump times
                pVal(3,t) = p;
                [h,p] = kstest2(testData3,jSizes,'Alpha',alpha); % test jump sizes
                pVal(4,t) = p;
            else
                pVal(3,t) = NaN;
                pVal(4,t) = NaN;
            end
        end
        
    case '3-OU-I1'
        
        for t = 1:Nmax/step
            t2 = t*step + burnin;      % index
            testData  = exprnd(par(t2,8),10000,1);
%             testData2 = exprnd( Tmax/par(t2,5) ,10000,1);
            %               compoundPoisson(Tmax,  t0,period,    delta,eta,jParameters,jDist)
            [~,testData2] = compoundPoisson(Tmax*5,par(t2,11),260,par(t2,10),par(t2,5),par(t2,8),2);
            testData2 = mod(testData2,260);
            testData3 = exprnd(par(t2,9),10000,1);
            testData4 = exprnd( Tmax/par(t2,7) ,10000,1);
            % ----------- Lchain1
            Phi = Lchain_1{t2};    % Phi iteration
            % get jump sizes
            jSizes = Phi(1,:);
            idx = jSizes>0;
            jSizes = jSizes(idx);
            % get jump times
            jTimes = Phi(2,idx);
            %jTimes = diff(jTimes);
            jTimes = mod(jTimes,260); % for inhomogenous Poisson process
            
            if numel(jTimes)>1                
                [h,p] = kstest2(testData2,jTimes,'Alpha',alpha); % test jump times                
                pVal(1,t) = p;
                [h,p] = kstest2(testData, jSizes,'Alpha',alpha); % test jump sizes
                pVal(2,t) = p;                                
            else
                pVal(1,t) = NaN;
                pVal(2,t) = NaN;
            end
            % ----------- Lchain2
            Phi = Lchain_2{t2};    % Phi iteration
            % get jump sizes
            jSizes = Phi(1,:);
            idx = jSizes>0;
            jSizes = jSizes(idx);
            % get jump times
            jTimes = Phi(2,idx);
            jTimes = diff(jTimes);
            if numel(jSizes)>1
                [h,p] = kstest2(testData4,jTimes,'Alpha',alpha); % test jump times
                pVal(3,t) = p;
                [h,p] = kstest2(testData3,jSizes,'Alpha',alpha); % test jump sizes
                pVal(4,t) = p;                
            else
                pVal(3,t) = NaN;
                pVal(4,t) = NaN;
            end
        end
end

pv = nanmean(pVal, 2); % compute mean p-value


end % end dJump3OU

