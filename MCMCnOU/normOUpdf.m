function pdf = normOUpdf(x, mu, lambda0, sigma)
%NORMOUPDF Compute the factors of the likelihood function of a Gaussian OU process
% 
% pdf = normOUpdf(x, mu, lambda0, sigma)
%
% The routine assumes that dt = 1
%
% Input arguments:
%
%   x       - vector of consecutive observations of the OU process
%
%   mu      - level of mean reversion
%
%   lambda0 - inverse of speed of mean reversion
%
%   sigma   - volatility
%
% Output arguments:
%
%   pdf     - vector containing the likelihood factors

dt = 1;     % diff(gridT);

% conditional mean of x_i given x_{i-1}
m = mu +(x(1:end-1) - mu).*exp(-dt/lambda0); 
% variance
v = sigma^2/2*lambda0*(1-exp(-2*dt/lambda0)) ; 
 % conditional pdf x_i given x_{i-1}
pdf = normpdf(x(2:end),m,sqrt(v));


