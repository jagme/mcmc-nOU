function [par, Lchain_, accRate] = gibbs(E, model, n, x0, hPar, propVar, NLupdates, str)
%GIBBS Hastings-withing-Gibbs algorithm for the 2-OU models
%
% [par, Lchain_, accRate] = gibbs(E, model, n, x0, hPar, propVar, NLupdates, str)
%
% Input arguments:
%
%   E  - vector of observations 
%
%   model - string to define model
%            '2-OU'     for 2-OU model: one positive jump component
%
%            '2-OU-I1'  for 2-OU-I1: one positive jump component with L a
%                       time inhomogeneous Poisson process
%
%            '2-OU-'    for 2-OU^{-}: one negative jump component                                   
%
%
%   n - number of MCMC samples to draw
%
%   x0 - structure array containing initial states for parameters
%
%   hPar - structure array containing prior hyperparameters
%
%   propVar - variance of proposal distributions
%
%   NLupdates - number of updates of process L (Phi) per MCMC iteration
%
%   str - string to save results, if empty no results are saved
%
% Output arguments
%
%   par             - Markov chain for all parameters
%
%   Lchain_         - Markov chain for marked Poisson process Phi
%
%   accRate         - number of accepted proposals for MH steps

Tmax = length(E) - 1;
period = 260;           % calendar year of 260 days

%------------------------------------------------
% proposal variances
sigmaPropLambda0 = propVar.lambda0;
sigmaPropLambda1 = propVar.lambda1;

%------------------------------------------------
% intial state jump process
L_dsizes = zeros(1,length(E)); %%% set L==0
Lchain0  = [L_dsizes; 0:Tmax];

%------------------------------------------------
% vectors for storing current and previous states
% will store output chain in different vectors
mu = zeros(2,1);        
sigma = zeros(2,1);
lambda0 = zeros(2,1);  lambda1 = zeros(2,1); 
gamma = zeros(2,1); % gamma = eta*Tmax, average number of jumps
beta  = zeros(2,1);

Lchain = cell(1,2);     % chain for processes L
% set initial states
mu(:) = x0.mu;
sigma(:) = x0.sigma;
lambda0(:) = x0.lambda0;   lambda1(:) = x0.lambda1;  
gamma(:) = x0.eta*Tmax;
beta(:) = x0.beta;

Lchain{1} = Lchain0;
Lchain{2} = Lchain0;

if strcmp(model,'2-OU-I1')
    theta1 = zeros(2,1);    theta2 = zeros(2,1);
    t0 = zeros(2,1);
    theta1(:) = x0.theta1;  theta2(:) = x0.theta2;
    t0(:) = x0.t00;
end

%------------------------------------------------
n2 = 20000;             % # of steps to save to show intermediate results
j = 0; j2 = 0;          % to manage intermediate results
before = 1; now = 2;    % index for the 2-size vectors
j3 = 100;               % store results every j3 steps
Lchain_ = cell(1,n/j3); % store L process
printR = 30000:50000:n; % print results at these steps + n2
%printR = 10000:20000:n;
k = 1; % index for printR

% run algorithm
switch model % using cases to optimise code for each model
    
    case '2-OU' % 2-OU model
        N_parameters = 6;
        parameters = zeros(n2,N_parameters);    % store parameters for intermediate results
        par = zeros(n/j3,N_parameters);         % store parameters every j3 steps                
        negJ = 1;                               % jump component Y1 is positive, w1 = 1
        accRate = zeros(5,1);                   % acceptance rates
        gibbs2OU;
        
    case '2-OU-' % 2-OU^{-} model
        N_parameters = 6;
        parameters = zeros(n2,N_parameters);    % store parameters for intermediate results
        par = zeros(n/j3,N_parameters);         % store parameters every j3 steps                
        negJ = -1;                              % jump component Y1 is negative, w1 = -1
        accRate = zeros(5,1);                   % acceptance rates
        gibbs2OU;
    
    case '2-OU-I1' % 2-OU-I1 model
        N_parameters = 6 + 2;
        parameters = zeros(n2,N_parameters);    % store parameters for intermediate results
        par = zeros(n/j3,N_parameters);         % store parameters every j3 steps
        negJ = 1;
        accRate = zeros(8,1);                   % acceptance rates
        gibbs2OUI1;
        
    otherwise
        error('Model not supported')
        
end

%------------------------------------------------ 
% save results
T = 0:Tmax;
Lchain = Lchain_;
% save results?
if ~isempty(str)
    save(str,'par','Lchain','T','E','accRate','hPar','x0', 'propVar', 'NLupdates','-v7.3')
end
    
% with Exponential jump sizes and natural parameterisation
    function gibbs2OU
        for i = 2:n
            if 0==mod(i,10000)
                fprintf('MCMC step %d out of %d\n',i,n);
            end
                        
            LT = Lchain{before};            
            XL = E - negJ*getY2(lambda1(before),LT,Tmax+1);
            
            %------------ update mu
            mu(now) = update_mu(XL, lambda0(before), sigma(before), hPar.Amu, hPar.Bmu);
            
            %------------ update sigma
            sigma(now) = update_sigma(XL, lambda0(before), mu(now), hPar.Asigma, hPar.Bsigma);
            
            %------------ update lambda0
            [lambda0(now), temp]= update_lambda0(XL,lambda0(before),mu(now),sigma(now), sigmaPropLambda0);
            accRate(1) = accRate(1) + temp;
                                             
            %------------ update lambda1
            [lambda1(now), temp ]= update_lambda1(E,LT,lambda1(before),lambda0(now),mu(now),sigma(now), sigmaPropLambda1,negJ);
            accRate(2) = accRate(2) + temp;     
            
            %--------------- update beta and gamma ( gamma = eta*Tmax)            
            beta(now)  = update_jumpBeta(LT,hPar.Bbeta);            
            gamma(now) = update_gamma(LT, hPar.Beta);
            
            %--------------- update L
            % prepare jump parameters [distribution parameter1 parameter2]
            % Exponential distribution
            jumpPar = [2 beta(now)];
            
            p = randi(2); % choose algorithm
            switch p
                case 1 % birth and death
                    for ii = 1:NLupdates
                    [Lchain(now), temp] = updateL_BD(E,LT,lambda0(now),mu(now),sigma(now),lambda1(now),jumpPar,gamma(now),negJ);                    
                    accRate(3) = accRate(3) + temp;
                    LT = Lchain{now};
                    end
                
                case 2 % displacement
                    [Lchain(now), temp] = updateL_dis(E,LT,lambda0(now),mu(now),sigma(now),lambda1(now),jumpPar,gamma(now),negJ);
                    accRate(4) = accRate(4) + temp;
            end
            
            %------------ update all jump sizes
            for ii = 1:NLupdates
            [Lchain(now), temp] = updateL_sizeBlocked(E,Lchain{now},lambda0(now),mu(now),sigma(now),lambda1(now),jumpPar,gamma(now),negJ);
            accRate(5) = accRate(5) + temp;
            end
            %NLT = length(LT(1,LT(1,:)>0)); % number of jumps           
            
            %------------ display intermediate results
            if ~isempty(printR)
                if (and(i>printR(k),i<=printR(k)+n2))
                    j = j+1;
                    %Njumps = sum(LT(1,LT(1,:)>0)); % sum of jump sizes                    
                    parameters(j,:) = [mu(now) lambda0(now) sigma(now) lambda1(now) gamma(now)/Tmax beta(now)];
                end
            end
            if j == n2
                j = 0; % reset counter and print intermediate results
                fprintf('%10s %10s %10s %10s %10s %10s \n', 'mu', 'lambda0', 'sigma', 'lambda1', 'eta','beta')                
                fprintf('%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f  \n', mean(parameters(:,:)))
                fprintf('%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f  \n', std(parameters(:,:)))
                %disp([mean(parameters(:,:)); std(parameters(:,:))])
                k = min(k+1,length(printR)); % get next k
                %plotR(parameters2(:,4),1,[3 3 1],'lambda_2')
            end
            % save elements for output chain
            if 0 == mod(i,j3)
                j2 = j2 + 1;
                %Njumps = sum(LT(1,LT(1,:)>0));
                par(j2,:) = [mu(now) lambda0(now) sigma(now) lambda1(now) gamma(now) beta(now)];
                Lchain_(j2) = Lchain(now);
                
            end
            % swap before and now states
            now = 1 + mod(i,2);
            before = 1 + mod(i+1,2);
        end        
        
    end % end gibbs Exponential

% with Exponential jump sizes and natural parameterisation, 2-OU-I1 model
    function gibbs2OUI1
        for i = 2:n
            if 0 == mod(i, 10000)
                fprintf('MCMC step %d out of %d\n',i,n);
            end
                        
            LT = Lchain{before};            
            XL = E - negJ*getY2(lambda1(before),LT,Tmax+1);
            
            %------------ update mu
            mu(now) = update_mu(XL, lambda0(before), sigma(before), hPar.Amu, hPar.Bmu);
            
            %------------ update sigma
            sigma(now) = update_sigma(XL, lambda0(before), mu(now), hPar.Asigma, hPar.Bsigma);
            
            %------------ update lambda1
            [lambda0(now), temp]= update_lambda0(XL,lambda0(before),mu(now),sigma(now), sigmaPropLambda0);
            accRate(1) = accRate(1) + temp;
                                             
            %------------ update lambda2
            [lambda1(now), temp ]= update_lambda1(E,LT,lambda1(before),lambda0(now),mu(now),sigma(now), sigmaPropLambda1,negJ);
            accRate(2) = accRate(2) + temp;     
            
            %--------------- update beta and gamma ( gamma = eta*Tmax)            
            beta(now) = update_jumpBeta(LT,hPar.Bbeta);            
%             gamma(now) = update_gamma(LT, hPar.Bgamma);
            % instead of updating gamma, we have three new parameters
            % (delta, eta, theta) = theta1, theta2, t0
            Tau = LT(2,LT(1,:)>0);
            [theta1(now),theta2(now),t0(now), temp] = updateIntensityFunction(Tmax,Tau,t0(before),period,theta1(before),theta2(before));
            accRate(6:8) = accRate(6:8) + temp;
            %--------------- update L
            % prepare jump parameters [distribution parameter1 parameter2]
            % Exponential distribution
            jumpPar = [2 beta(now)];
            gamma_now = struct('eta',0,'t0',t0(now),'period',period,'theta1',theta1(now),'theta2',theta2(now));
            p = randi(2); % choose algorithm
            switch p
                case 1 % birth and death
                    for ii=1:NLupdates
                    [Lchain(now), temp] = updateL_BD_2OUI1(E,LT,lambda0(now),mu(now),sigma(now),lambda1(now),jumpPar,gamma_now, negJ);                    
                    accRate(3) = accRate(3) + temp;
                    LT = Lchain{now};
                    end
                
                case 2 % displacement
                    [Lchain(now), temp] = updateL_dis2OUI1(E,LT,lambda0(now),mu(now),sigma(now),lambda1(now),jumpPar,gamma_now, negJ);
                    accRate(4) = accRate(4) + temp;
            end
            % update all jump sizes with RW
            for ii = 1:NLupdates
            [Lchain(now),temp] = updateL_sizeBlocked(E,Lchain{now},lambda0(now),mu(now),sigma(now),lambda1(now),jumpPar,NaN,negJ);
            accRate(5) = accRate(5) + temp;
            end
            %NLT = length(LT(1,LT(1,:)>0));
            
            %intermediate results
            if ~isempty(printR)
                if (and(i>printR(k),i<=printR(k)+n2))
                    j = j+1;
                    %Njumps = sum(LT(1,LT(1,:)>0));
                    parameters(j,:) = [mu(now) lambda0(now) sigma(now) lambda1(now) beta(now) theta1(now) theta2(now) t0(now)];
                end
            end
            if j == n2
                j = 0; % reset counter and print intermediate results
                fprintf('%10s %10s %10s %10s %10s %10s %10s %10s \n', 'mu', 'lambda0', 'sigma', 'lambda1', 'beta','delta1', 'eta1', 'theta1')
                fprintf('%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n', mean(parameters(:,:)))
                fprintf('%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n', std(parameters(:,:)))
                %disp([mean(parameters(:,:));std(parameters(:,:))])
                k = min(k+1,length(printR)); % get next k
                %plotR(parameters2(:,4),1,[3 3 1],'lambda_2')
            end
            % save elements for output chain
            if 0 == mod(i,j3)
                j2 = j2+1;
                %Njumps = sum(LT(1,LT(1,:)>0));
                par(j2,:) = [mu(now) lambda0(now) sigma(now) lambda1(now) beta(now) theta1(now) theta2(now) t0(now)];
                Lchain_(j2) = Lchain(now);
                
            end
            % swap before and now states
            now = 1 + mod(i,2);
            before = 1 + mod(i+1,2);
        end
        
    end % end gibbs Exponential 2-OU-I1


end % end main function