function [newL, k] = updateL_BD_3OUI1(X, oldL, lambda1, mu, sigma, lambda2, jumpPar, eta, Laux, lambdaaux, negJ)
%UPDATEL_BD_3OUI1 Birth-and-death algorithm for the marked Poisson process Phi_i of the 3-OU-I1 model 
%
% [newL, k] = updateL_BD_3OUI1(X, oldL, lambda1, mu, sigma, lambda2, jumpPar, eta, Laux, lambdaaux, negJ)
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
%   jumpPar - current state of parameters of jump distribution
%           jumpPar(1) is the distribution
%           jumpPar(2:end) is the parameters
%
%   eta - structure array containing current states of jump intensity function 
%           parameters: eta, delta, theta, and the fixed period
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
T = length(X) - 1;
rate = eta.eta/T;
t0 = eta.t0;
period = eta.period;
theta1 = eta.theta1;
theta2 = eta.theta2;
jParameter = jumpPar(2:end);
jumpDist = jumpPar(1);

newL = cell(1,1);

k = 0;      % # accepted proposals

p = 0.5;    % probability of birth move

L1 = oldL;
idx_locations = find(L1(1,:)>0);
N_locations = length(idx_locations);

r = 0;
if rand <= p % generate new point
    c = rand*T; % new-born time
    jump = drawJumpSize(jParameter,jumpDist,1); % new jump
    L2 = [L1 [jump;c]];
    L2 = sortrows(L2',2)';% sort L2 wrt to time variable
    
    if negJ(1) == 1 % if updating first jump component
        rate = intensityFun(c,t0,period,theta1,theta2);
    end
    r = likelihood*rate*T*(1-p)/(p*(N_locations+1));
    
elseif N_locations>0
    % remove one point from process
    a = randi(N_locations,1); % choose random point
    L2 = L1;
    L2(:,idx_locations(a)) = []; % remove it
    c = L1(2,idx_locations(a));
    if negJ(1) == 1 % if updating first jump component
        rate = intensityFun(c,t0,period,theta1,theta2);
    end
    temp = rate*T*(1-p)/(p*(N_locations+1));
    r = likelihood/temp;
end

% accept or reject proposal
if rand <= r % accept
    newL{1} = L2;
    k = k + 1;
else
    newL{1} = L1; % keep initial L
end

    function l = likelihood
        % likelihood ratio        
        Y22 = getY2(lambdaaux,Laux,lenX);
        XL1 = X - negJ(1)*getY2(lambda2,L1,lenX) - negJ(2)*Y22;
        XL2 = X - negJ(1)*getY2(lambda2,L2,lenX) - negJ(2)*Y22;
        
        L_L = prod(normOUpdf(XL2,mu,lambda1,sigma)./normOUpdf(XL1,mu,lambda1,sigma));
        
        l = L_L;
        
    end
end