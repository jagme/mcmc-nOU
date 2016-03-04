function plotPhi(file1, file2)
%PLOTPHI Plot true process Phi = (Phi1, Phi2) from a path of the 3-OU model
%and the corresponding MCMC-estimated process
%
% Input arguments:
%
%   file1 - file path to true processes LT1 and LT2 obtained from function
%       simulate3OUmodel.m
%
%   file2 - file path to MCMC results for the estimated 3-OU model
%

load(file1)

load(file2)

%%
% set time range to plot
Lint = 1;                          % left index
Rint = Lint + length(E) - 1;       % right index, the whole range in this case
Nk = min(5000, length(Lchain_1));  % number of mcmc Phi elements to use, Nk < # length of chain
T_ = T(Lint:Rint);
% Tmax = max(T_);
% posterior summary for L1
Lchain = Lchain_1;
F1 = summaryL(Lchain(end-Nk+1:end),T_);  % get estimate of Phi1
% repeat for L2
Lchain = Lchain_2;
F2 = summaryL(Lchain(end-Nk+1:end),T_);  % get estimate of Phi2

% process estimates
meanL1  = summaryL2(F1,T_,Nk,0.4);
meanL2  = summaryL2(F2,T_,Nk,0.4);

%% plot results
figure;
row = 4; col= 1;
xlimY = 4;
% % ------------------ L1
subplot(row,col,1)
bar(LT1(2,:),LT1(1,:),'BarWidth',0.4,'FaceColor','red','EdgeColor','blue')
xlim([T_(1) T_(end)]), ylim([0 xlimY])
xlabel('Time')
title('True process \Phi_1')
grid on
% ------------------
subplot(row,col,2)
bar(meanL1(2,:),meanL1(1,:),'BarWidth',.4,'FaceColor','b','EdgeColor','blue'), 
xlim([T_(1) T_(end)]), ylim([0 xlimY])
title('Estimated \Phi_1')
xlabel('Time')
hold on
grid on

% % ------------------ L2
subplot(row,col,3)
bar(LT2(2,:),LT2(1,:),'BarWidth',.4,'FaceColor','red','EdgeColor','blue') 
xlim([T_(1) T_(end)]),ylim([0 xlimY])
xlabel('Time')
title('True process \Phi_2')
grid on

subplot(row,col,4)
bar(meanL2(2,:),meanL2(1,:),'BarWidth',.4,'FaceColor','b','EdgeColor','blue'), 
xlim([T_(1) T_(end)]), ylim([0 xlimY])
xlabel('Time')
hold on
grid on
title('Estimated \Phi_2')