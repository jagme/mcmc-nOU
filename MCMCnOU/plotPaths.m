function plotPaths(file1,file2)
%PLOTPATHS Plot paths of the estimated 2-OU and 3-OU models 
%for the APXUK 2001 - 2006 data as shown in Figure 3
%
% Input arguments:
%
%   file1 - path to MCMC results for the 2-OU model
%
%   file2 - path to MCMC results for the 3-OU model

% 2-OU 

load(file1)         % load file for 2-OU model

% set time range
Lint = 451;
Rint = Lint + 200;
T_ = T(Lint:Rint);
% Tmax = max(T_);
% Ltrial = Lchain{end};

figure;
subplot(2,1,1)                  
% plot deseasonalised data
plot(T_,E(Lint:Rint))
hold on
% get Y1
Y22 = getY2(par(end,4),Lchain{end},length(E));
% get Y0
Y0 = E - Y22;
% plot Y0
plot(T_,Y0(Lint:Rint))

grid on, box off
legend({'Deseasonalised APXUK','Path of Y_0'},'FontSize',11, 'Location','northwest'),
xlabel('Time (days)'), ylabel('£/MWh')
title('Estimated paths of the 2-OU model for the APXUK index')
h = gca; h.FontSize = 11;

subplot(2,1,2)
% plot Y1
plot(T_, Y22(Lint:Rint),'Color',[0.93 0.69 0.13])
grid on, box off
xlim([T_(1) T_(end)])
legend({'Path of Y_1'},'FontSize',11,'Location','northwest')
xlabel('Time (days)'), ylabel('£/MWh')
h = gca; h.FontSize = 11;

%% 3-OU
load(file2)             % load file for 3-OU

t2 = length(par);
Ltrial1 = Lchain_1{end};
Ltrial2 = Lchain_2{end};
Y1 = getY2(par(t2,4),Ltrial1,length(E));
Y2 = getY2(par(t2,6),Ltrial2,length(E));
Y0 = E - Y1 - Y2;

T_ = [450 650];

figure;

subplot(2,1,1)
% plot deseasonalised data
plot(E) 
title('Estimated paths of the 3-OU model for the APXUK index')
hold on, 
% plot Y0
plot(Y0)
legend({'Deseasonalised APXUK','Path of Y_0'},'FontSize',11,'Location','northwest')
xlabel('Time (days)'), ylabel('£/MWh'), xlim(T_)
box off, grid on
h = gca; h.FontSize = 11;

subplot(2,1,2)
% plot Y1 & Y2
plot(Y1, 'Color',[0.93 0.69 0.13])
xlim(T_)
hold on
plot(Y2, 'Color',[0.49 0.18 0.56])
legend({'Path of Y_1','Path of Y_2'},'FontSize',11,'Location','northwest')
xlim(T_), ylim([0 6])
box off, grid on
h = gca; h.FontSize = 11;
xlabel('Time (days)'), ylabel('£/MWh')