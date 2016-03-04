function [newL, k] = updateL_BD(X, oldL, lambda0, mu, sigma, lambda1, jumpPar, eta, negJ)
%UPDATEL_BD Birth-and-death algorithm for the marked Poisson process of the 2-OU model 
%
% [newL, k] = updateL_BD(X, prop, lambda0, mu, sigma, lambda1, jumpPar, eta, negJ)
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
%   lambda1 - current state of parameter lambda1
%
%   jumpPar - current state of parameters of jump distribution
%           jumpPar(1) is the distribution
%           jumpPar(2:end) is the parameters
%
%   eta - current state of gamma = eta * Tmax
%
%   negJ - indicates whether jump component Y1 is positive (negJ = 1) or
%   negative (negJ = -1)
%
% Output arguments:
%
%   newL - new state of process L
%
%   k - if proposal was accepted k = 1, otherwise k = 0
%
lenX = length(X);
T = length(X) - 1;
rate = eta/T;
jParameter = jumpPar(2:end);
jumpDist = jumpPar(1);

newL = cell(1,1);

k = 0;   % # accepted proposals
p = 0.5; % probability of birth move

L1 = oldL;
idx_locations = find(L1(1,:)>0);
N_locations = length(idx_locations);

r = 0;
if rand <= p    % generate new point
    c = rand*T; % new-born time
    jump = drawJumpSize(jParameter,jumpDist,1); % new jump
    L2 = [L1 [jump;c]];
    L2 = sortrows(L2',2)';% sort L2 wrt to time variable
    % likelihood ratio
    r = likelihood*rate*T*(1-p)/(p*(N_locations+1));
    
elseif N_locations>0
    % remove one point from process
    a = randi(N_locations,1); % choose random point
    L2 = L1;
    %L2(1,idx_locations(a)) = 0;
    L2(:,idx_locations(a)) = []; % remove it
    %c = L1(2,idx_locations(a));
    temp = rate*T*(1-p)/(p*(N_locations+1));
    r = likelihood/temp;
end

% accept or reject proposal
if rand <= r        % accept
    newL{1} = L2;
    k = k + 1;
else
    newL{1} = L1;  % keep initial L
end

    function l = likelihood
        % likelihood for birth-death move
        XL1 = X - negJ*getY2(lambda1,L1,lenX);
        XL2 = X - negJ*getY2(lambda1,L2,lenX);
        L_L = prod(normOUpdf(XL2,mu,lambda0,sigma)./normOUpdf(XL1,mu,lambda0,sigma));
        
        l = L_L;
        
    end
end