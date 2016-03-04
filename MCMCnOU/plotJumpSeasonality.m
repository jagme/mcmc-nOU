function plotJumpSeasonality(filepath, burnin)
%PLOTJUMPSEASONALITY Plot average number of jumps per year versus 
%theoretical jump intensity function I1 for the EEX 2000 - 2006, 
%as shown in Figure 3
%
% Input arguments:
%
% filepath - file path to 3-OU-I1 mcmc results from gibbs3OU
%
% burnin - burn in period for mcmc chain
%
% This file assumes that 
% EEX 2000 - 2006 dataset starts on June 16, 2000
% (APXUK starts on March 27, 2001)

load(filepath)

idxMonths = linspace(0,260,13); % time in months, calendar year of 260 days
Nmax = length(Lchain_1);
if burnin>=Nmax
    error('Burn-in period should be less than %d', Nmax)
end
idx2 = NaN(12,Nmax-burnin-1);
idx3 = NaN(12,Nmax-burnin-1);
for j = 1:Nmax-burnin-1
    Ltrial = Lchain_1{burnin+j};  % Phi_1 for iteration burn + j
    J = Ltrial(1,:);            
    idx = J>0.;                 % index with all jumps
    % the jump times
    % add phase to times to start on 1 January
    phase = 11 + 130; % EEX: 11 +130 , APXUK: 70 + 130
    T = mod(Ltrial(2,idx) + phase,260); 
    [jT,I] = sort(T);           % sort times
    J = J(idx); J = J(I);
    for i = 1:12
        temp = and(jT>=idxMonths(i),jT<idxMonths(i+1)); % jumps on each month
        idx2(i,j) = sum(temp);      % how many jumps in month i
        idx3(i,j) = mean(J(temp));  % mean jump size on month i
    end
end

% plot results

figure
bar(mean(idx2,2)/length(E)*260,'FaceColor',[0.3 0.75 0.93]) % yearly average
xlim([0 13])
% title('Average number of positive spikes on the EEX market'), %xlabel('Month')
box off, %ylim([0 45])
h = gca; h.FontSize = 11;
h.XTick =  1:12;
h.XTickLabel = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
h.XTickLabelRotation = 20;

%% plot theoretical intensity function
% parameters in yearly units 
theta  = 1 + mod(mean(par(burnin:end,11)),130)/260; % estimate of theta1 + 1, period of I1 is 130 days
period = 12;    % 12 months == 260 days 
delta  = mean(par(burnin:end,10));                  % estimate of delta1
eta    = mean(par(burnin:end,5));                   % estimate of eta1

Tmax = 12; % 12 months

t = linspace(1,Tmax,10000);
y = intensityFun(t,theta,period,delta,eta) * 5/.25;

hold on
ph = plot(t,y,'k');
xlim([0, 13])
% legend(ph,'Fitted intensity function \times 20')
legend({'Average','Fitted intensity function \times 20'},'FontSize',11)
