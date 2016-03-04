function [x1, k] = update_lambda2_3OU(E, LT2, x0, lambda0, mu, sigma, sigmaProp, lambda1, LT1, negJ)
%UPDATE_LAMBDA2_3OU MH step for inverse of speed of mean reversion lambda1 of 
% Non-Gaussian OU process Y1 in the 3-OU model
%
% [x1, k] = update_lambda2_3OU(E, LT2, x0, lambda0, mu, sigma, sigmaProp, lambda1, LT1, negJ)
%
% Input arguments:
%
%   E - observed data 
%
%   LT1 - current state of Poisson process L1
%
%   x0 - current state of lambda1
%
%   lambda0 - current state of lambda0
%
%   mu - current state of mu
%
%   sigma - current state of sigma
%
%   sigmaProp - SD for proposal distribution
%
%   lambda2 - current state of inverse of speed of mean reversion of Y2
%
%   LT2 - current state of Poisson process L2
%
%   negJ - vector of size 2 indicating sign of jump components
%
% Output arguments:
%
%   x1 - new state of lambda1
%
%   k - logical value to indicate whether the proposed value was accepted 
%       1 if proposal was accepted, 0 otherwiseMCMC for relaxation time lambda2 in Y2
% multiplicative random walk MH algorithm

Tmax = length(E) - 1;

x = exp(-1/x0);  
y = exp(normrnd(log(x),sigmaProp));

proposal = -1/log(y);

Y1L1 = getY2(lambda1,LT1,Tmax+1);
XL_x = E - negJ(1)*Y1L1 - negJ(2)*getY2(x0,LT2,Tmax+1);         % get observations of Y1
XL_y = E - negJ(1)*Y1L1 - negJ(2)*getY2(proposal,LT2,Tmax+1);   % get observations of Y1

% if both jump components have the same sign impose the restriction 
% lambda1 > lambda2
if negJ(1) == negJ(2)
    if proposal < lambda1
        L_L = prod(normOUpdf(XL_y,mu,lambda0,sigma)./normOUpdf(XL_x,mu,lambda0,sigma))...
            * betapdf(y,1,1)/betapdf(x,1,1) *y/x;
    else
        L_L = 0;
    end
else
    L_L = prod(normOUpdf(XL_y,mu,lambda0,sigma)./normOUpdf(XL_x,mu,lambda0,sigma))...
            * betapdf(y,1,1)/betapdf(x,1,1) *y/x;
end
if rand <= L_L
    x1 = proposal;
    k = 1;
else
    x1 = x0;
    k = 0;
end


