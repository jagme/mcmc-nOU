function [newL, k] = updateL_sizeBlocked_3OU(X, oldL, lambda0, mu, sigma, lambdai, jumpPar, Laux, lambdaaux, negJ)
%UPDATEL_SIZEBLOCKED_3OU Block update of jump sizes for the marked Poisson process
%Phi_i of the 3-OU models
%
% [newL, k] = updateL_sizeBlocked_3OU(X, oldL, lambda0, mu, sigma, lambda1, jumpPar, Laux, lambdaaux, negJ)
%
% Input arguments:
%
%   X - observed data
%
%   oldL - current state of process Phi_i
%
%   lambda0, mu, sigma - current states of parameters lambda0, mu & sigma
%   in Y0
%
%   lambdai - current state of parameter lambdai
%
%   jumpPar - current state of parameters of jump distribution for oldL
%           jumpPar(1) is the distribution
%           jumpPar(2:end) is the parameters
%
%   Laux - the marked Poisson process of the other jump component, if
%   updating Phi_1, oldL = Phi_1 and Laux = Phi_2 and vice versa
% 
%   lambdaaux - the inverse of speed of mean reversion corresponding to the
%   OU process having Laux
%   
%   negJ - indicates the sign of the jump components, the first element
%   correspond to the component being simulated
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
    Yaux = getY2(lambdaaux,Laux,lenX);
    XL1 = X - negJ(1)*getY2(lambdai,L1,lenX) - negJ(2)*Yaux;
    XL2 = X - negJ(1)*getY2(lambdai,L2,lenX) - negJ(2)*Yaux;
    % if exponential jump size
%     switch jumpDist
%         case 1 % Pareto
            %temp = (oldjump/newjump)^jParameter(2);
%         case 2 % exponential
            temp = prod(factorJump)*exp(-1/jParameter*sum(newjump-oldjump));
%     end
    % ratio
    L_L = prod(normOUpdf(XL2,mu,lambda0,sigma)./normOUpdf(XL1,mu,lambda0,sigma))...
        * temp;
    r = L_L;
end

% accept or reject proposal
if rand <= r % accept
    newL{1} = L2;
    k = k + 1;
else
    newL{1} = L1; % keep initial L
end

end