function x1 = update_sigma(Z, lambda0, mu, a, b)
%UPDATE_SIGMA Gibbs sampler for sigma
%
% x1 = update_sigma(Z, lambda0, mu, a, b)
%
% Input arguments:
%
%   Z       - deflated observations Z = X - \sum Yi
%
%   lambda0 - current state of inverse of speed of mean reversion of Y0
%
%   mu      - current state of mu
%
%   a, b    - hyperparameters for Inverse-Gamma(a,b) prior
%
% Output arguments:
%
%   x1 - sample of sigma
%
n = length(Z);
dt = 1;
temp  = exp(-1/lambda0*dt);
temp2 = mu*(1-temp);
A = n/2;
B = 1/lambda0/(1-temp^2)*sum((Z(2:end)-temp2-Z(1:end-1)*temp).^2);

x1 = gamrnd(A+a,1/(B+b),1,1); % gamma(A+a,B+b) mean A+b / B+b

x1 = sqrt(1/x1); % return sqrt of Inverse-Gamma sample, i.e sample of sigma