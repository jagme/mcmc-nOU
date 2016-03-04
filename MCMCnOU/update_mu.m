function x1 = update_mu(Z, lambda0, sigma, Amu, Bmu)
%UPDATE_MU Gibbs sampler for level of mean reversion mu of Y0
%
% x1 = update_mu(Z, lambda0, sigma, Amu, Bmu)
%
% Input arguments
%
%   Z       - deflated observations Z = X - \sum Yi
%
%   lambda0 - current state of inverse of speed of mean reversion of Y0
%
%   sigma   - current state of sigma
%
%   Amu, Bmu   - hyperparameters for Normal prior N(Amu, Bmu^2)
%

n = length(Z);

speedMR = 1/lambda0;    % speed of mean reversion
dt = 1;                 % assumed timestep
temp = exp(-speedMR*dt);
v2 = lambda0*sigma^2*(1-temp^2)/2; % variance OU process
temp2 = 1 - temp;
temp3 = v2/Bmu^2;
% flat prior
%A = sum((XL(2:end)-XL(1:end-1)*exp(-speedMR*dt)))/n/(1-temp);
%B = sigma^2*(1-temp^2)/2/speedMR/n/(1-temp)^2;
% normal prior
A = (temp2*sum((Z(2:end)-Z(1:end-1)*exp(-speedMR*dt))) + Amu*temp3) / (n*temp2^2 +temp3 );
B = v2/(temp2^2*n+temp3);

x1 = normrnd(A, sqrt(B), 1, 1);
