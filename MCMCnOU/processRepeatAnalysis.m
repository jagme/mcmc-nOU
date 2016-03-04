function ci = processRepeatAnalysis(strpath, burnin)
%PROCESSREPEATANALYSIS Process the results of function repeatAnalysis 
% and print averages of posterior means alongside their 'confidence intervals'
%
% Input arguments:
%
%   strpath - file path to saved results
%
%   burnin - burn-in period for chain
%

load(strpath) % load results
SizeChain = size(chain);
if SizeChain(3) <= burnin
   error('Burn-in period should be less than %d', SizeChain(3)) 
end

R  = SizeChain(2); % number of mcmc chains generated
Me = SizeChain(1); % number of different parameters eta
Np = SizeChain(4); % number of parameters 
est = NaN(4,R,Np);
est2 = NaN(Np+2,Np+2);           % mean of posterior means
ci = NaN(Np+2,Me*3);             % 'confidence intervals'

for i = 1:Me
    for j = 1:R    
        est(i,j,:) = mean(chain(i,j,burnin:end,:));
    end
    temp2 = squeeze(est(i,:,:));
    temp = temp2;
    temp(:,2) = temp2(:,3).^2;               % sigma^2
    temp(:,[3 5]) = exp(-1./temp2(:,[2 4])); % rho0, rho1
    temp(:,[4 6]) = temp2(:,[2 4]);          % lambda0, lambda1
    temp(:,7) = temp2(:,5) / 1000;           % eta
    temp(:,8) = temp2(:,6);                  % beta
    
    est2(:,(i-1)*2+1:i*2) = [mean(temp)' std(temp)'];    
    ci(:,(i-1)*3+1:i*3) = [mean(temp)' mean(temp)'-1.96*std(temp)' mean(temp)'+1.96*std(temp)'];
%     ci(:,(i-1)*3+1:i*3) = [mean(temp)' mean(temp)'-1.96*std(temp)'/sqrt(R) mean(temp)'+1.96*std(temp)'/sqrt(R)];

    
end

strpar = {'mu','sigma^2', 'rho0', 'lambda0', 'rho1', 'lambda1','eta','beta'};
fprintf('Average and spread of posterior means\n')
fprintf('----------------------------------------------------------------------------------------------\n')
fprintf('True eta      eta = %.3f  \t eta = %.3f   \t      eta = %.3f \t    eta = %.3f\n', eta_)
fprintf('----------------------------------------------------------------------------------------------\n')
for i = 1:Np+2    
        fprintf('%10s %10.4f %20.4f %20.4f %20.4f\n', strpar{i}, ci(i, 1:3:end))
        fprintf('%17.4f, %1.4f %12.4f, %1.4f %12.4f, %1.4f %12.4f, %1.4f\n', ci(i,[2 3 5 6 8 9 11 12]))    
end



