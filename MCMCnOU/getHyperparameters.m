function hPar = getHyperparameters(model)
% Convenience function to set prior hyperparameters for each model
%
% hPar = getHyperparameters(model)
%
switch model
    case {'2-OU','2-OU-' '2-OU-I1'}
        hPar.Amu = 1;           % mu
        hPar.Bmu = 20;
        
        hPar.Asigma = 1.5;      % sigma
        hPar.Bsigma = 0.005;        
        
        hPar.Bbeta = 1;         % beta
        
        hPar.Beta = 10;         % eta
                
    case {'3-OU', '3-OU-', '3-OU-I1'}
        hPar.Amu = 1;           % mu
        hPar.Bmu = 20;
        
        hPar.Asigma = 1.5;      % sigma
        hPar.Bsigma = 0.005;
        
        hPar.Bbeta1 = 1;        % beta1
        
        hPar.Bbeta2 = 1;        % beta2
        
        hPar.Beta1 = 10;        % eta1
        
        hPar.Beta2 = 10;        % eta2
        
end
