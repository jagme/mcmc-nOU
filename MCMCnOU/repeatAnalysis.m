function chain = repeatAnalysis(filepath)
%REPEATANALYSIS Repeated estimation of the 2-OU model. 
% This routine calls the MCMC algorithm N times with different sample paths 
% of the the 2-OU model and saves each output chain. With function 
% processRepeatAnalysis.m compute mean of posterior means and confidence
% intervals.
%
% Input arguments:
%
%   filepath - file path name to save results
%

T    = 1000;                % Tmax for simulation of 2-OU model
eta_ = [0.05 0.1 0.2 0.3 ]; % true jump rate, eta
Me   = length(eta_);          
N    = 60;                  % number of calls to gibbs
% burn = 1;                 % burn-in period
% parM = NaN(Me,N,6);
n = 1000000;                % number of MCMC iterations

chain = NaN(Me,N,n/100,6);

for i = 1:Me
    eta = eta_(i);
    parfor j = 1:N                
        par = changePar(T,eta,n);
%         parM(i,j,:) = mean(par(burn:end,:));
        chain(i,j,:,:) = par;
    end
    fprintf('Finished i = %d out of %d', i, Me)
end

% 'results/2c-estimatesSimT1000.mat'
save(filepath, 'chain','eta_','T')
%mean(parM)