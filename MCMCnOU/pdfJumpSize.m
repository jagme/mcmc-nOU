function pdf = pdfJumpSize(x,parameters,distribution)
%PDFJUMPSIZE Convenience function to compute pdf of a jump size
%
% pdf = pdfJumpSize(x,parameters,distribution)
%
% Input arguments:
%
%   x   - jump size
%
%   parameters - parameters of jump size distribution
%
%   distribution - the jump size distribution
%               1: Pareto(z0, alpha)
%               2: Exponential(lambda)
%
% Output arguments:
%
% pdf - pdf of jump size using specified probability distribution
%

switch distribution
    % Pareto (z0,alpha)
    case 1        
        z0 = parameters(1); alpha = parameters(2);
        pdf = alpha*z0^alpha./x.^(alpha+1);
    % exponential (lambda)        
    case 2        
        pdf = exppdf(x,parameters);
end