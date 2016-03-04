function [x1, k] = update_lambda1(E, LT, x0, lambda0, mu, sigma, sigmaProp, negJ)
%UPDATE_LAMBDA1 MH step for inverse of speed of mean reversion lambda1 of 
% Non-Gaussian OU process Y1 in the 2-OU model
%
% [x1, k] = update_lambda1(E, LT, x0, lambda0, mu, sigma, sigmaProp, negJ)
%
% Input arguments:
%
%   E - observed data 
%
%   LT - current state of Poisson process L
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
%   negJ - indicate whether jump component is positive (negJ = 1) or
%   negative (negJ = -1)
%
% Output arguments:
%
%   x1 - new state of lambda1
%
%   k - logical value to indicate whether the proposed value was accepted 
%       1 if proposal was accepted, 0 otherwise


Tmax = length(E) - 1;

x = exp(-1/x0);  y = exp(normrnd(log(x),sigmaProp));
proposal = -1/log(y);

XL_x = E - negJ*getY2(x0,LT,Tmax+1);
XL_y = E - negJ*getY2(proposal,LT,Tmax+1);

L_L = prod(normOUpdf(XL_y,mu,lambda0,sigma)./normOUpdf(XL_x,mu,lambda0,sigma))...
    * betapdf(y,1,1)/betapdf(x,1,1) *y/x;

if rand < L_L
    x1 = proposal;
    k = 1;
else
    x1 = x0;
    k = 0;
end


