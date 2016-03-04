function plotR(chain,burn,place,str_par)
%PLOTR Plot summary of MCMC chain: trace, ACF and posterior density
%
% plotR(chain,burn,place,str_par)
%
% Input arguments:
%
% chain     - chain of parameter
%
% burn      - burn-in period, the first 'burn' samples are removed for
%               computations
%
% place     - if using subplot specify the row for this panel of graphs
%
% str_par   - string with name of parameter
%

% remove first elements of chain
chain = chain(burn:end);
% plot chain
subplot(place(1),place(2),place(3))
plot(chain),ylabel(str_par), xlim([0,length(chain)])
hold on,
% plot mean
plot(1:500,mean(chain),'r')
% plot posterior density
subplot(place(1),place(2),place(3)+1)
ksdensity(chain), xlabel(str_par), %ylabel('Density'),
%xlim([min(chain),max(chain)])
grid on
% plot ACF
nLags = min(1001, length(chain)-1);
try
    [ACF,lags,bounds] = autocorr(chain,nLags,nLags-1,[]);  
    % idx = find(ACF<bounds(1));
    % if ~isempty(idx) idx = round(idx(1)/50+1)*50+100;else idx = nLags; end
    idx = nLags;
    lags = lags(1:idx); ACF = ACF(1:idx) ;
    h = subplot(place(1),place(2),place(3)+2);
    plot(lags,ACF)
    xlim([0,lags(end)])
    % ylim([bounds(2)*2, 1])
    ylim([bounds(2) - 0.1, 1])
    hold on
    line([lags(1) lags(end)],[bounds(1) bounds(1)],'Color','k','LineStyle','--')
    line([lags(1) lags(end)],[bounds(2) bounds(2)],'Color','k','LineStyle','--')
    %ylabel('ACF'), xlabel('Lag'),
    grid on
    %set(h,'XTick',0:.2:1)
    %set(h,'XTickLabel',{'1'})
catch Me
    switch Me.identifier
        case 'MATLAB:UndefinedFunction'
            warning(Me.message)
        otherwise
        rethrow(Me)
    end
end