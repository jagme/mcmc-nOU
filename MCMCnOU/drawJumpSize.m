function jumpsize = drawJumpSize(parameters,distribution,n)
%DRAWJUMPSIZE draw jump sizes variates
%
% jumpsize = drawJumpSize(parameters,distribution,n)
%
% Input arguments:
%
%   parameters - parameters for jump size distribution
%
%   distribution -  a number specifying the distribution
%    1: Pareto
%    2: Exponential with mean parameters(1)
%    3: Gamma with mean parameters(1)/parameter(2) 

if nargin<3
    n = 1;
end

switch distribution
    
    % Pareto (z0,alpha)
    case 1        
        z0 = parameters(1); 
        alpha = parameters(2);
        jumpsize = z0./rand(1,n).^(1/alpha);

    % Exponential (lambda)
    case 2  % exp with mean parameter parameters(1)              
        jumpsize = exprnd(parameters(1),1,n);
    
    % Gamma(a,b) with mean a/b
    case 3
        jumpsize = gamrnd(parameters(1),1/parameters(2));
end
