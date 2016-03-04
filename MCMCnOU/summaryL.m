function F = summaryL(chain,T)
%SUMMARYL Summarise posterior/Markov chain of compound Process L (Phi)
%
% F = summaryL(chain,T)
%
% Input arguments:
%
%   chain - cell array containing chain for process L, ie this is obtained
%   as output of MCMC algorithm: gibbs or gibbs3OU
%
%   T - vector of time points defining the time intervals to summarise L
%       
% Output arguments:
%
%   F - matrix containing summary of the posterior of L
%   For each interval i
%   F(i, :) = [interval number-of-jumps std-of-jump-sizes number-unique-jumps mean-jump-time]

F = NaN(length(T)-1,6); % contains summary

for i = 1:length(T)-1 % for each interval
%     if mod(i,50)==0
%         disp(i)
%     end
    leftend = T(i); rightend = T(i+1); % end points of interval
    temp = [];
    temp2 = [];
    for j = 1:length(chain) % for each MCMC sample
        
        aux_L = chain{j};
        % find jumps in the half-open interval (leftend,rightend]
        idx1 = find(and(aux_L(2,:)>leftend  ,aux_L(2,:)<=rightend));
        idx2 = aux_L(1,idx1)>0;         % find jumps
        idx3 = idx1(idx2);              % indexes of jumps on interval
        temp = [temp aux_L(1,idx3)];    % accumulate jump sizes
        temp2 = [temp2 aux_L(2,idx3)];  % accumulate jump times
        
    end % end loop over elements
    
    if isempty(temp)
        F(i,:) = [i 0 0 0 0 T(i)];
    else % interval, number of jumps, mean, std, # unique jumps, mean time
        F(i,:) = [i length(temp) mean(temp) std(temp) length(unique(temp)) mean(temp2)];
    end
end

F = [0 0 0 0 0 0;F]; % F(i) contains info about interval i = 1,2,...