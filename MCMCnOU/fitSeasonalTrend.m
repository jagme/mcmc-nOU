function [E, beta] = fitSeasonalTrend(datastr, plotFit)
%FITSEASONALTREND Fit seasonal trend function to dataset specified by datastr
%
% [E, beta] = fitSeasonalTrend(datastr, plotFit)
%
% Input arguments:
%
%   datastr  - string specifying dataset
%               'apxuk1'  APXUK 27/03/2001 - 21/11/2006
%               'apxuk2'  APXUK 24/01/2011 - 16/02/2015
%               'eex1'    EEX   16/06/2000 - 21/11/2006 
%               'eex2'    EEX   24/01/2011 - 16/02/2015
%
%   plotFit - boolean variable, if true results are plotted
%
%   Output arguments:
%
%   E  - deseasonalised series
%
%   beta - fitted parameters
%%
switch datastr  
    case 'apxuk1'        
        idx1 = 1; idx2 = 1476;         % 27/03/2001 - 21/11/2006
    case 'apxuk2'        
        idx1 = 2565; idx2 = 3625;      % 24/01/2011 - 16/02/2015
    case 'eex1'        
        idx1 = 1; idx2 = 1678;         % 16/06/2000 - 21/11/2006
    case 'eex2'        
        idx1 = 2767; idx2 = 3827;      % 24/01/2011 - 16/02/2015
    otherwise
        error('datastr should be either apxuk1, apxuk2, eex1 or eex2')
end

load('../data/APXUK-EEXdatasets.mat');

if strfind(datastr,'apx')
    namedata = 'APXUK';
    units = '£';
    E0 = APXPowerUKSpotBaseLoadIndexTimeSeriesData.Price;
else
    namedata = 'EEX';
    units = '?';
    E0 = EEXBASETimeSeriesDataeuro.Price;
end

E0 = E0(idx1:idx2)';
E1 = E0;

% there are 3 large negative prices on EEX2, average them out with neighbours 
idx = find(E0<0);
if ~isempty(idx)
    E0(idx(1:2)) = mean(E0([idx(1)-1 idx(2)+1]));
    E0(idx(3))   = mean(E0([idx(3)-1 idx(3)+1]));
end

% specify seasonal trend function
period = 260; % period 365 days - weekends
strendFun = @(a,t)(a(1) + a(2)*t +a(3)*sin(2*pi*t/period)+ a(4)*cos(2*pi*t/period)...
    + a(5)*sin(4*pi*t/period)+a(6)*cos(4*pi*t/period));

% initial guess
beta0 = [1 1 1 1 1 1];
t = 1:length(E0);
% set robust options for nlinfit
% options = statset('nlinfit');
% options.Robust = 'on';
% fit
[beta, R]= nlinfit(t,log(E0),strendFun,beta0);
% deseasonalise raw data
E = E1./exp(strendFun(beta,t));

% plot results?
if plotFit    
    
    figure;    
    plot(log(E0))    
    hold on 
    grid on
    plot(strendFun(beta,t),'black')            
    plot(E), grid on, box off
    
    xlabel('Time (days)'), ylabel([units '/MWh'])
    legend(['Daily mean ' namedata ' log-prices'],'Seasonal trend', ['Deseasonalised ' namedata ' series' ])    
    h = gca; h.FontSize = 11;
end

