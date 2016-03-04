function [chain, k] = update_lambda0(XL,x0, mu,sigma, sigmaProp)
%UPDATE_LAMBDA0 MH step for inverse of speed of mean reversion lambda0
%
% [chain, k] = update_lambda1(XL,x0, mu,sigma, sigmaProp)
%
% Input arguments:
%
%   XL - deflated data Z = X - Y1
%
%   x0 - current state of lambda0
%
%   mu - current state of mu
%
%   sigma - current state of sigma
%
%   sigmaProp - SD for proposal distribution
%
% Output arguments:
%
%   x1 - new state of lambda0
%
%   k - logical value to indicate whether the proposed value was accepted 
%       1 if proposal was accepted, 0 otherwise

x = exp(-1/x0);
y = exp(normrnd(log(x),sigmaProp));
proposal = -1/log(y);
% MH ratio
L_L = prod(normOUpdf(XL,mu,proposal,sigma)./normOUpdf(XL,mu,x0,sigma))...
    * betapdf(y,1,1)/betapdf(x,1,1) * y/x;

if rand < L_L
    chain = proposal;
    k = 1;
else
    chain = x0;
    k = 0;
end