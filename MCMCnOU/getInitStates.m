function x0 = getInitStates(model)
% Convenience function to set initial state of chain for each model
%
% x0 = getInitStates(model)
%

model = upper(model);

switch model
    
    case {'2-OU','2-OU-'}
        
        x0.mu = 1;
        
        x0.lambda0 = 5;
        
        x0.sigma = 0.1;
        
        x0.lambda1 = 2;
        
        x0.eta = 0.001;
        
        x0.beta = 0.5;                                
    
    case '2-OU-I1'
        
        x0.mu = 1;
        
        x0.lambda0 = 5;
        
        x0.sigma = 0.1;
        
        x0.lambda1 = 2;
        
        x0.eta = 0.001;
        
        x0.beta = 0.5;
        
        x0.theta1 = 0.5; % delta
        
        x0.theta2 = 0.5; % eta
        
        x0.t00 = 100;    % theta
                
    case {'3-OU','3-OU-'}
        
        x0.mu = 1;
        
        x0.lambda0 = 5;
        
        x0.sigma = 0.2;
        
        x0.lambda1 = 5;
        x0.lambda2 = 1;
        
        x0.eta1 = 0.001;
        x0.eta2 = 0.001;
        
        x0.beta1 = 0.5;
        x0.beta2 = 0.5;                
        
    case '3-OU-I1'
        
        x0.mu = 1;
        
        x0.lambda0 = 5;
        
        x0.sigma = 0.2;
        
        x0.lambda1 = 5;
        x0.lambda2 = 1;
        
        x0.eta1 = 0.001;
        x0.eta2 = 0.001;
        
        x0.beta1 = 0.5;
        x0.beta2 = 0.5;
        
        x0.theta1 = 0.5; % delta
        
        x0.theta2 = 0.5; % eta
        
        x0.t00 = 100;    % theta
        
end