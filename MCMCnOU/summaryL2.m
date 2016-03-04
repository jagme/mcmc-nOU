function meanL = summaryL2(F,T_,Nk,prob)
%SUMMARYL2 Remove jumps whose frequency is too small in Markov chain
%
% meanL = summaryL2(F,T_,Nk,prob)
%
% Input arguments:
%
%   F   - output from summaryL.m
%
%   T_  - vector of time points defining the time intervals to summarise L
%
%   Nk  - number of mcmc samples used to get F
%
%   prob - remove jumps whos relative frequency across MCMC iterations is
%       less than prob
%
% Output arguments:
%
% meanL - 'posterior summary' of L
%

temp = F(:,2)/Nk;   %max(F(:,2)); % # jumps in interval / max # jumps in all intervals
idx = temp<prob;    % plot mean jump size whose frequency is greater than prob
F(idx,2:end-1) = 0; % erase jumps with small frequency
F(idx,end) = floor(F(idx,end)); % floor of mean jump times

% get a proper process L from the summary F, based on mean of jump sizes and mean of
% jump times for each interval

%         mean size,  mean time
meanL = [F(2:end,3) F(2:end,6) ]';  %zeros(2,length(T_));
missL = setdiff(T_,F(2:end,6));     % get missing times
meanL = [meanL [zeros(1,length(missL));missL ]];
meanL = sortrows(meanL',2)';
meanL(2,:) = meanL(2,:) - T_(1); % need time from 0
