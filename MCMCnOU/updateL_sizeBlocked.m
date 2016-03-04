function [newL, k] = updateL_sizeBlocked(X, oldL, lambda0, mu, sigma, lambda1, jumpPar, eta, negJ)
%UPDATEL_sizeBlocked Block update of jump sizes of marked Poisson process
%Phi of the 2-OU model
%
% [newL, k] = updateL_sizeBlocked(X, oldL, lambda0, mu, sigma, lambda1, jumpPar, eta, negJ)
%
% Input arguments:
%
%   X - observed data
%
%   oldL - current state of process L
%
%   lambda0, mu, sigma - current states of parameters lambda0, mu & sigma
%   in Y0
%
%   lambda1 - current state of parameter lambda1 associated with Y1
%
%   jumpPar - current state of parameters of jump distribution
%           jumpPar(1) is the distribution
%           jumpPar(2:end) is the parameters
%
%   eta - current state of gamma = eta * Tmax
%   
%   negJ - indicates the sign of the jump component Y1
%
% Output arguments:
%
%   newL - new state of process L
%
%   k - if proposal was accepted k = 1, otherwise k = 0
%

lenX = length(X);
%T = length(X)-1;
%rate = eta/T;
jParameter = jumpPar(2:end);
jumpDist = jumpPar(1);

newL = cell(1,1);

k = 0; % # accepted proposals
% generate chain

L1 = oldL;
idx_locations = find(L1(1,:)>0);
N_locations = length(idx_locations);

r = 0;
if N_locations>0
    oldjump = L1(1,idx_locations); % all jumps
    
    temp = 0.8; % to control variance of multiplicative Random Walk    
    variance = 2*log(1-log(temp)/N_locations/4);    
    factorJump = exp(randn(1,N_locations)*sqrt(variance)); % one for each jump
    newjump = oldjump.*factorJump;
    L2 = L1;
    L2(1,idx_locations) = newjump; % change jump
     
    XL1 = X - negJ*getY2(lambda1,L1,lenX);
    XL2 = X - negJ*getY2(lambda1,L2,lenX);
    
    % if exponential jump size
%     switch jumpDist
%         case 1 % Pareto
            %temp = (oldjump/newjump)^jParameter(2);
%         case 2 % exponential
            temp = prod(factorJump)*exp(-1/jParameter*sum(newjump-oldjump));
%     end
    % likelihood ratio
    L_L = prod(normOUpdf(XL2,mu,lambda0,sigma)./normOUpdf(XL1,mu,lambda0,sigma))...
        * temp;
    r = L_L;
end

% accept or reject proposal
if rand <= r        % accept proposal
    newL{1} = L2;
    k = k + 1;
else
    newL{1} = L1;  % keep initial L
end

end
