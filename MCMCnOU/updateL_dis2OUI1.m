function [newL, k] = updateL_dis2OUI1(X, oldL, lambda0, mu, sigma, lambda1, jumpPar, sPar, negJ)
%UPDATEL_DIS2OUI1 Local displacement move for the marked Poisson process of
%2-OU-I1 model
%
% [newL, k] = updateL_BD_2OUI1(X, oldL, lambda1, mu, sigma, lambda2, jumpPar, eta)
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
%           jumpPar(2:end) is/are the parameter(s)
%
%   eta - structure array with current states intensity function
%   parameters: eta1, delta1, theta1 and fixed period
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
%rate = eta/T;
t0 = sPar.t0;
period = sPar.period;
theta1 = sPar.theta1;
theta2 = sPar.theta2;
jParameter = jumpPar(2:end);
jumpDist = jumpPar(1);

newL = cell(1,1);

k = 0;  % # accepted proposals
% generate chain

L1 = oldL;
idx_locations = find(L1(1,:)>0);
N_locations = length(idx_locations);
r = 0;

if N_locations>0
    a = randi(N_locations,1); % choose random point/time
    location =idx_locations(a);
    if a==1 % need to be careful with endpoints
        tau_i = 0;
    else
        tau_i = L1(2,idx_locations(a-1));
    end
    if a==N_locations
        tau_i2=T;
    else
        tau_i2 = L1(2,idx_locations(a+1));
    end
    newtau = tau_i +  rand*(tau_i2-tau_i); % new time is unif[tau(i-1), tau(i+1)]
    oldtau = L1(2,location);
    oldjump = L1(1,location);
    C = 1;
    newjump = oldjump*exp(-1/lambda1*(newtau - oldtau)/C);
%     if jumpDist==1, newjump=max(jParameter(1),newjump); end
    L2 = L1;
    L2(:,idx_locations(a)) = [newjump;newtau]; % displace point
    L2 = sortrows(L2',2)';% sort L2 wrt to time variable
    
    % likelihood and acceptance probability
    XL1 = X - negJ*getY2(lambda1,L1,lenX);
    XL2 = X - negJ*getY2(lambda1,L2,lenX);
    
%     switch jumpDist
%         case 1 % Pareto
%             temp = (oldjump/newjump)^jParameter(2) * exp(-1/lambda1*(newtau-oldtau)/C);
%         case 2 % exponential
            %temp = exp(jParameter*(oldjump-newjump) - 1/lambda2*(newtau-oldtau)/C);
            temp = exp(jParameter*(oldjump-newjump) - 1/lambda1*(newtau-oldtau)/C) ...
                * intensityFun(newtau,t0,period,theta1,theta2)/intensityFun(oldtau,t0,period,theta1,theta2);
%     end
    
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