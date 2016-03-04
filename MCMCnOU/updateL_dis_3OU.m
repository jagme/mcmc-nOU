function [newL, k] = updateL_dis_3OU(X, oldL, lambda0, mu, sigma, lambdai, jumpPar, Laux, lambdaaux, negJ)
%UPDATEL_DIS_3OU Local displacement move for the marked Poisson process Phi_i of the 3-OU model 
%
% [newL, k] = updateL_dis_3OU(X, oldL, lambda0, mu, sigma, lambda1, jumpPar, Laux, lambdaaux, negJ)
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
%   lambdai - current state of parameter lambdai associated with Yi
%
%   jumpPar - current state of parameters of jump distribution
%           jumpPar(1) is the distribution
%           jumpPar(2:end) is the parameters
%
%   eta - current state of gammai = etai * Tmax
%
%   Laux - the marked Poisson process of the other jump component, if
%   updating Phi_1, oldL = Phi_1 and Laux = Phi_2 and vice versa
% 
%   lambdaaux - the inverse of speed of mean reversion corresponding to the
%   OU process having Laux as Levy process
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
    C=1;
    newjump = oldjump*exp(-1/lambdai*(newtau - oldtau)/C);
    if jumpDist==1, newjump=max(jParameter(1),newjump); end
    L2 = L1;
    L2(:,idx_locations(a)) = [newjump;newtau]; % displace point
    L2 = sortrows(L2',2)';% sort L2 wrt to time variable
        
    Y22 = getY2(lambdaaux,Laux,lenX);
    XL1 = X - negJ(1)*getY2(lambdai,L1,lenX) - negJ(2)*Y22;
    XL2 = X - negJ(1)*getY2(lambdai,L2,lenX) - negJ(2)*Y22;
    
%     switch jumpDist
%         case 1 % Pareto
%             temp = (oldjump/newjump)^jParameter(2) * exp(-1/lambdai*(newtau-oldtau)/C);
%         case 2 % exponential
            temp = exp(jParameter*(oldjump-newjump) - 1/lambdai*(newtau-oldtau)/C);
%     end
    L_L = prod(normOUpdf(XL2,mu,lambda0,sigma)./normOUpdf(XL1,mu,lambda0,sigma))...
        * temp;
    r = L_L;
end

% accept or reject proposal
if rand <= r        % accept
    newL{1} = L2;
    k = k+1;
else
    newL{1} = L1;  % keep initial L
end

end