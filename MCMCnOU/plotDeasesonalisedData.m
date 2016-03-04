function plotDeasesonalisedData
% Plot Figure 1

% get deseasonalised data
E_APXUK1 = fitSeasonalTrend('apxuk1', false);
E_APXUK2 = fitSeasonalTrend('apxuk2', false);
E_EEX1   = fitSeasonalTrend('eex1',   false);
E_EEX2   = fitSeasonalTrend('eex2',   false);

%%
Nfill = 250;

figure;

% panel 1
subplot(2, 1, 1)
Euk = [E_APXUK1 Inf(1, Nfill) E_APXUK2];    % the APXUK data
plot(Euk)
legend('Deseasonalised APXUK'), ylabel('£/MWh') % \pounds / MWh
grid on, box off
h = gca; h.FontSize = 11;
h.XTick = [ 12 142:130:(1678+29) (1678+29+Nfill-24):130:3000];
h.XTickLabel = {'','Jan 2001','','Jan 2002','','Jan 2003',...
    '','Jan 2004','','Jan 2005','','Jan 2006','','Jan 2007','Jan 2011','',...
    'Jan 2012','','Jan 2013','','Jan 2014','','Jan 2015'};
h.XTickLabelRotation = 20;

% panel 2
subplot(2, 1, 2)
Eeex = [E_EEX1 Inf(1, Nfill) E_EEX2];
plot(Eeex)
legend('Deseasonalised EEX'), ylabel('?/MWh') % \euro / MWh
grid on, box off
h = gca; h.FontSize = 11;
h.XTick = [ 12 142:130:(1678+29) (1678+29+Nfill-24):130:3000];
h.XTickLabel = {'','Jan 2001','','Jan 2002','','Jan 2003',...
    '','Jan 2004','','Jan 2005','','Jan 2006','','Jan 2007','Jan 2011','',...
    'Jan 2012','','Jan 2013','','Jan 2014','','Jan 2015'};
h.XTickLabelRotation = 20;
