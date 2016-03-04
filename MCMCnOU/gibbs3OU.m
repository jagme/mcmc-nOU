function [par, accRate] = gibbs3OU(E, model, n, x0, hPar, propVar, NLupdates, str)
%GIBBS3OU Hastings-within-Gibbs algorithm for the 3-OU models
%
% [par, accRate] = gibbs3OU(E, model, n, x0, hPar, propVar, N_Lupdates, str)
%
% Input arguments:
%
%   E  - vector of observations
%
%   model - string to define model
%            '3-OU'     for 3-OU model: two positive jump components
%
%            '3-OU-I1'  for 3-OU-I1: one positive jump component with L1 a
%                       time inhomogeneous Poisson process, one negative
%                       jump component. w1 = 1, w2 = -1
%
%            '3-OU-'    for 3-OU^{-}: one positive and one negative jump 
%                       component. w1 = 1, w2 = -1
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
%   accRate         - number of accepted proposals for MH steps

Tmax = length(E) - 1;
period = 260;

%------------------------------------------------
% proposal variances
sigmaPropLambda0 = propVar.lambda0;
sigmaPropLambda1 = propVar.lambda1;
sigmaPropLambda2 = propVar.lambda2;

%------------------------------------------------
% intial states jump process
L_dsizes = zeros(1,length(E)); %%% set L==0
Lchain0  = [L_dsizes;0:Tmax];
Lchain30 = [L_dsizes;0:Tmax];

%------------------------------------------------
% vectors for storing current and previous states
% will store output chain in different vectors
mu    = zeros(2,1);       
sigma = zeros(2,1);
lambda0 = zeros(2,1); lambda1 = zeros(2,1); lambda2 = zeros(2,1);
gamma1  = zeros(2,1); gamma2  = zeros(2,1); % gammai = etai*Tmax
beta1   = zeros(2,1); beta2   = zeros(2,1);
% for processes Phi
Lchain1 = cell(1,2);   Lchain2 = cell(1,2);
% set initial states
mu(:)    = x0.mu; 
sigma(:) = x0.sigma;
lambda0(:) = x0.lambda0;    lambda1(:) = x0.lambda1; lambda2(:) = x0.lambda2; 
gamma1(:)   = x0.eta1 * Tmax; gamma2(:)   = x0.eta2 * Tmax;
beta1(:)   = x0.beta1;       beta2(:)   = x0.beta2;
Lchain1{1} = Lchain0;        Lchain1{2}  = Lchain0;
Lchain2{1} = Lchain30;       Lchain2{2}  = Lchain30;

if strfind(model,'3-OU-I1')
    theta1 = zeros(2,1);    theta2 = zeros(2,1);
    t0 = zeros(2,1);
    theta1(:) = x0.theta1;  theta2(:) = x0.theta2;
    t0(:) = x0.t00;
end

%------------------------------------------------
n2 = 20000;          % # of steps to save to show intermediate results

j = 0; j2 = 0;       % to manage intermediate results
before = 1; now = 2; % index for the 2-size vectors
% save results every j3 steps
j3 = 100;
Lchain_1 = cell(1,n/j3); % and L process here
Lchain_2 = cell(1,n/j3); % and L process here
printR = 30000:50000:n;  % print results at these steps + n2
%printR = 10000:20000:n;
k = 1; % index for printR

% run algorithm
switch model % using this to optimise code for each model    
    case '3-OU' % Exponential
        N_parameters = 9;
        parameters = zeros(n2,N_parameters);    % store parameters for intermediate results
        par = zeros(n/j3,N_parameters); % store parameters every 100 steps        
        negJ = [1 1]; % both jump components are positive, w1 = 1; w2 = 1;        
        accRate = zeros(9,1);   % acceptance rate for MH steps
        gibbs3OU;
    case '3-OU-' % 3-OU-n, Y1 positive Y2 negative
        N_parameters = 9;
        parameters = zeros(n2,N_parameters);    % store parameters for intermediate results
        par = zeros(n/j3,N_parameters); % store parameters every 100 steps
        negJ = [1 -1]; % both jump components are positive, w1 = 1; w2 = -1;
        accRate = zeros(9,1);   % acceptance rate for MH steps
        gibbs3OU;
    case '3-OU-I1' % 3-OU-I1, Y1 positive Y2 negative, L2 time-inhomogeneous
        N_parameters = 11;
        parameters = zeros(n2,N_parameters);    % store parameters for intermediate results
        par = zeros(n/j3,N_parameters); % store parameters every 100 steps
        negJ = [1 -1]; % both jump components are positive, w1 = 1; w2 = -1;
        accRate = zeros(12,1);   % acceptance rate for MH steps
        gibbs3OUI1;
    % cases for prior sensitivity
    case '3-OU-I1-sensitivity-sigma'
        N_parameters = 11;
        parameters = zeros(n2,N_parameters);    % store parameters for intermediate results
        par = zeros(n/j3,N_parameters); % store parameters every 100 steps
        negJ = [1 -1]; % both jump components are positive, w1 = 1; w2 = -1;
        accRate = zeros(13,1);   % acceptance rate for MH steps
        gibbs3OUI1Sigma;
        
    case '3-OU-I1-sensitivity-eta'
        N_parameters = 11;
        parameters = zeros(n2,N_parameters);    % store parameters for intermediate results
        par = zeros(n/j3,N_parameters); % store parameters every 100 steps
        negJ = [1 -1]; % both jump components are positive, w1 = 1; w2 = -1;
        accRate = zeros(13,1);   % acceptance rate for MH steps
        gibbs3OUI1Eta;
        
        
end


%------------------ save results
T = 0:Tmax;

if ~isempty(str)
    save(str,'par','Lchain_1','Lchain_2','T','E','hPar','x0','propVar','accRate','NLupdates','-v7.3')
end

% MCMC algorithm for 3-OU models
    function gibbs3OU
        for i = 2:n
            if 0==mod(i,10000)
                fprintf('MCMC step %d out of %d\n',i,n);
            end
            
            LT1 = Lchain1{before};
            LT2 = Lchain2{before};
            XL = E - negJ(1)*getY2(lambda1(before),LT1,Tmax+1) - negJ(2)*getY2(lambda2(before),LT2,Tmax+1);                        
            
            %------------ update mu
            mu(now) = update_mu(XL, lambda0(before), sigma(before), hPar.Amu, hPar.Bmu);
            
            %------------ update sigma
            sigma(now) = update_sigma(XL, lambda0(before), mu(now), hPar.Asigma, hPar.Bsigma);
            
            %------------ update lambda0
            [lambda0(now), temp]= update_lambda0(XL,lambda0(before),mu(now),sigma(now), sigmaPropLambda0);
            accRate(1) = accRate(1) + temp;
            
            %------------ update lambda1
            [lambda1(now), temp] = update_lambda1_3OU(E,LT1,lambda1(before),...
                lambda0(now),mu(now),sigma(now), sigmaPropLambda1, lambda2(before), LT2,negJ);
            accRate(2) = accRate(2) + temp;
                        
            %------------ update lambda2
            [lambda2(now), temp] = update_lambda2_3OU(E,LT2,lambda2(before),...
                lambda0(now),mu(now),sigma(now), sigmaPropLambda2, lambda1(now), LT1,negJ);
            accRate(3) = accRate(3) + temp;
                       
            %------------ update beta1
            beta1(now) = update_jumpBeta(LT1,hPar.Bbeta1);
            
            %------------ update beta2
            beta2(now) = update_jumpBeta(LT2,hPar.Bbeta2);
            
            %------------ update eta1
            gamma1(now) = update_gamma(LT1,hPar.Beta1);
            
            %------------ update eta2
            gamma2(now) = update_gamma(LT2,hPar.Beta2);
            
            %------------ update L1
            % prepare jump parameters [distribution parameter1 parameter2]           
            jumpPar = [2 beta1(now)];                        
            p = randi(2);
            switch p
                case 1 % birth and death  
                    for ii = 1:NLupdates
                    [Lchain1(now), temp] = updateL_BD_3OU(E,LT1,lambda0(now),mu(now),sigma(now),lambda1(now),jumpPar,gamma1(now),LT2,lambda2(now),negJ);                    
                    LT1 = Lchain1{now};
                    accRate(4) = accRate(4) + temp;
                    end
                case 2 % displacement
                    [Lchain1(now), temp] = updateL_dis_3OU(E,LT1,lambda0(now),mu(now),sigma(now),lambda1(now),jumpPar,LT2,lambda2(now),negJ);
                    accRate(5) = accRate(5) + temp;
            end
            % update all jump sizes with RW
            for ii = 1:NLupdates
            [Lchain1(now),temp] = updateL_sizeBlocked_3OU(E,Lchain1{now},lambda0(now),mu(now),sigma(now),lambda1(now),jumpPar,LT2,lambda2(now),negJ);
            accRate(6) = accRate(6) + temp;
            end
            %------------ update L2
            jumpPar = [2 beta2(now)];
            p = randi(2);%p=0;
            negJflip = [negJ(2) negJ(1)];
            switch p
                case 1 % birth and death  
                    for ii = 1:NLupdates
                    [Lchain2(now), temp] = updateL_BD_3OU(E,LT2,lambda0(now),mu(now),sigma(now),lambda2(now),jumpPar,gamma2(now),LT1,lambda1(now),negJflip);                    
                    LT2 = Lchain2{now};
                    accRate(7) = accRate(7) + temp;
                    end
                case 2 % displacement
                    [Lchain2(now),temp] = updateL_dis_3OU(E,LT2,lambda0(now),mu(now),sigma(now),lambda2(now),jumpPar,LT1,lambda1(now),negJflip);
                    accRate(8) = accRate(8) + temp;
            end
            for ii = 1:NLupdates
            Lchain2(now) = updateL_sizeBlocked_3OU(E,Lchain2{now},lambda0(now),mu(now),sigma(now),lambda2(now),jumpPar,Lchain1{now},lambda1(now),negJflip);
            accRate(9) = accRate(9) + temp;
            end            
            
            %------------ Display intermediate results
            if ~isempty(printR)
                if (and(i>printR(k),i<=printR(k)+n2))
                    j = j+1;                    
                    parameters(j,:) = [mu(now) lambda0(now) sigma(now) lambda1(now) gamma1(now)/Tmax lambda2(now) gamma2(now)/Tmax beta1(now) beta2(now)];
                end
            end
            if j == n2
                j = 0; % reset counter and print intermediate results
                fprintf('%10s %10s %10s %10s %10s %10s %10s %10s %10s \n', ...
                    'mu', 'lambda0', 'sigma', 'lambda1', 'eta1','lambda2', 'eta2','beta1', 'beta2')
                fprintf('%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', mean(parameters(:,:)))
                fprintf('%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', std(parameters(:,:)))
                %disp([mean(parameters(:,:));std(parameters(:,:))])
                k = min(k+1,length(printR)); % get next k
                %plotR(parameters2(:,4),1,[3 3 1],'lambda_2')
            end
            % save elements for output chain
            if 0 == mod(i,j3)
                j2 = j2 + 1;                
                par(j2,:) = [mu(now) lambda0(now) sigma(now) lambda1(now) gamma1(now) lambda2(now) gamma2(now) beta1(now) beta2(now) ];
                Lchain_1(j2) = Lchain1(now);
                Lchain_2(j2) = Lchain2(now);
            end
            % swap before and now states
            now = 1 + mod(i,2);
            before = 1 + mod(i+1,2);
        end
                
    end % end gibbs Exponential

% Algorithm for 3-OU-I1
    function gibbs3OUI1
        for i = 2:n
            if 0==mod(i,10000)
                fprintf('MCMC step %d out of %d\n',i,n);
            end
            
            LT1 = Lchain1{before};
            LT2 = Lchain2{before};
            XL = E - negJ(1)*getY2(lambda1(before),LT1,Tmax+1) - negJ(2)*getY2(lambda2(before),LT2,Tmax+1);                        
            
            %------------ update mu
            mu(now) = update_mu(XL, lambda0(before), sigma(before), hPar.Amu, hPar.Bmu);
            
            %------------ update sigma
            sigma(now) = update_sigma(XL, lambda0(before), mu(now), hPar.Asigma, hPar.Bsigma);
            
            %------------ update lambda0
            [lambda0(now), temp]= update_lambda0(XL,lambda0(before),mu(now),sigma(now), sigmaPropLambda0);
            accRate(1) = accRate(1) + temp;
            
            %------------ update lambda1
            [lambda1(now), temp] = update_lambda1_3OU(E,LT1,lambda1(before),...
                lambda0(now),mu(now),sigma(now), sigmaPropLambda1, lambda2(before), LT2,negJ);
            accRate(2) = accRate(2) +temp;
                        
            %------------ update lambda2
            [lambda2(now), temp] = update_lambda2_3OU(E,LT2,lambda2(before),...
                lambda0(now),mu(now),sigma(now), sigmaPropLambda2, lambda1(now), LT1,negJ);
            accRate(3) = accRate(3) + temp;
                       
            %------------ update beta1
            beta1(now) = update_jumpBeta(LT1,hPar.Bbeta1);
            
            %------------ update beta2
            beta2(now) = update_jumpBeta(LT2,hPar.Bbeta2);
            
            %------------ update eta1
            % gamma(now) = update_gamma(LT1,hPar.Beta1);
            
            %------------ update eta2
            gamma2(now) = update_gamma(LT2,hPar.Beta2);
            
            %------------ update intensity function parameters
            Tau = LT1(2,LT1(1,:)>0);
            [theta1(now),theta2(now),t0(now),temp] = updateIntensityFunction(Tmax,Tau,t0(before),period,theta1(before),theta2(before));
            accRate(10:12) = accRate(10:12) + temp;
            %------------ update L1
            % prepare jump parameters [distribution parameter1 parameter2]           
            jumpPar = [2 beta1(now)];    
            gamma_now = struct('eta',0,'t0',t0(now),'period',period,'theta1',theta1(now),'theta2',theta2(now));
            p = randi(2);
            switch p
                case 1 % birth and death  
                    for ii = 1:NLupdates
                    [Lchain1(now), temp] = updateL_BD_3OUI1(E,LT1,lambda0(now),mu(now),sigma(now),lambda1(now),jumpPar,gamma_now,LT2,lambda2(now),negJ);                    
                    LT1 = Lchain1{now};
                    accRate(4) = accRate(4) + temp;
                    end
                case 2 % displacement
                    [Lchain1(now), temp] = updateL_dis_3OUI1(E,LT1,lambda0(now),mu(now),sigma(now),lambda1(now),jumpPar,LT2,lambda2(now),negJ,gamma_now);
                    accRate(5) = accRate(5) + temp;
            end
            % update all jump sizes with RW
            for ii = 1:NLupdates
            [Lchain1(now), temp] = updateL_sizeBlocked_3OU(E,Lchain1{now},lambda0(now),mu(now),sigma(now),lambda1(now),jumpPar,LT2,lambda2(now),negJ);
            accRate(6) = accRate(6) + temp;
            end
            %------------ update L2
            jumpPar = [2 beta2(now)];
            gamma3_now = struct('eta',gamma2(now),'t0',t0(now),'period',period,'theta1',theta1(now),'theta2',theta2(now));
            p = randi(2);
            negJflip = [negJ(2) negJ(1)];
            switch p
                case 1 % birth and death  
                    for ii = 1:NLupdates
                    [Lchain2(now), temp] = updateL_BD_3OUI1(E,LT2,lambda0(now),mu(now),sigma(now),lambda2(now),jumpPar,gamma3_now,LT1,lambda1(now),negJflip);                    
                    LT2 = Lchain2{now};
                    accRate(7) = accRate(7) + temp;
                    end
                case 2 % displacement
                    [Lchain2(now),temp] = updateL_dis_3OUI1(E,LT2,lambda0(now),mu(now),sigma(now),lambda2(now),jumpPar,LT1,lambda1(now),negJflip,gamma3_now);
                    accRate(8) = accRate(8) + temp;
            end
            for ii = 1:NLupdates
            [Lchain2(now),temp] = updateL_sizeBlocked_3OU(E,Lchain2{now},lambda0(now),mu(now),sigma(now),lambda2(now),jumpPar,Lchain1{now},lambda1(now),negJflip);
            accRate(9) = accRate(9) + temp;
            end
            
            %intermediate results
            if ~isempty(printR)
                if (and(i>printR(k),i<=printR(k)+n2))
                    j = j+1;                    
                    parameters(j,:) = [mu(now) lambda0(now) sigma(now) lambda1(now) theta2(now) lambda2(now) gamma2(now)/Tmax beta1(now) beta2(now) theta1(now) t0(now)];
                end
            end
            if j == n2
                j = 0; % reset counter and print intermediate results
                fprintf('%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n', ...
                    'mu', 'lambda0', 'sigma', 'lambda1','eta1','lambda2', 'eta2','beta1', 'beta2','delta1','theta1')
                fprintf('%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', mean(parameters(:,:)))
                fprintf('%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', std(parameters(:,:)))
%                 disp([mean(parameters(:,:));std(parameters(:,:))])
                k = min(k+1,length(printR)); % get next k
                %plotR(parameters2(:,4),1,[3 3 1],'lambda_2')
            end
            % save elements for output chain
            if 0 == mod(i,j3)
                j2 = j2 + 1;                
                par(j2,:) = [mu(now) lambda0(now) sigma(now) lambda1(now) theta2(now) lambda2(now) gamma2(now) beta1(now) beta2(now) theta1(now) t0(now)];
                Lchain_1(j2) = Lchain1(now);
                Lchain_2(j2) = Lchain2(now);
            end
            % swap before and now states
            now = 1 + mod(i,2);
            before = 1 + mod(i+1,2);
        end
                
    end % end gibbs3OUI1

% Algorithm for 3-OU-I1 with flat prior for sigma^2
    function gibbs3OUI1Sigma
        for i = 2:n
            if 0==mod(i,10000)
                fprintf('MCMC step %d out of %d\n',i,n);
            end
            
            LT1 = Lchain1{before};
            LT2 = Lchain2{before};
            XL = E - negJ(1)*getY2(lambda1(before),LT1,Tmax+1) - negJ(2)*getY2(lambda2(before),LT2,Tmax+1);                        
            
            %------------ update mu
            mu(now) = update_mu(XL, lambda0(before), sigma(before), hPar.Amu, hPar.Bmu);
            
            %------------ update sigma
            % prior for sigma is Uniform(0, 0.25)
%             sigma(now) = update_sigma(XL, lambda_1(before), mu(now), hPar.Asigma, hPar.Bsigma);
            x = sigma(before);  y = exp(normrnd(log(x),0.05));                                       
            L_L = prod(normOUpdf(XL,mu(now),lambda0(now),y)./normOUpdf(XL,mu(now),lambda0(now),x))...
                     * y/x * (y<0.25)*(x<0.25);                                    
            if rand <= L_L
               sigma(now) = y;
               accRate(13) = accRate(13) + 1;
            else
               sigma(now) = x; 
            end
            %------------ update lambda0
            [lambda0(now), temp]= update_lambda0(XL,lambda0(before),mu(now),sigma(now), sigmaPropLambda0);
            accRate(1) = accRate(1) + temp;
            
            %------------ update lambda1
            [lambda1(now), temp] = update_lambda1_3OU(E,LT1,lambda1(before),...
                lambda0(now),mu(now),sigma(now), sigmaPropLambda1, lambda2(before), LT2,negJ);
            accRate(2) = accRate(2) + temp;
                        
            %------------ update lambda2
            [lambda2(now), temp] = update_lambda2_3OU(E,LT2,lambda2(before),...
                lambda0(now),mu(now),sigma(now), sigmaPropLambda2, lambda1(now), LT1,negJ);
            accRate(3) = accRate(3) + temp;
                       
            %------------ update beta1
            beta1(now) = update_jumpBeta(LT1,hPar.Bbeta1);
            
            %------------ update beta2
            beta2(now) = update_jumpBeta(LT2,hPar.Bbeta2);
            
            %------------ update eta1
            % gamma(now) = update_gamma(LT1,hPar.Beta1);
            
            %------------ update eta2
            gamma2(now) = update_gamma(LT2,hPar.Beta2);
            
            %------------ update intensity function parameters
            Tau = LT1(2,LT1(1,:)>0);
            [theta1(now),theta2(now),t0(now),temp] = updateIntensityFunction(Tmax,Tau,t0(before),period,theta1(before),theta2(before));
            accRate(10:12) = accRate(10:12) + temp;
            
            %------------ update L1
            % prepare jump parameters [distribution parameter1 parameter2]           
            jumpPar = [2 beta1(now)];    
            gamma_now = struct('eta',0,'t0',t0(now),'period',period,'theta1',theta1(now),'theta2',theta2(now));
            p = randi(2);
            switch p
                case 1 % birth and death  
                    for ii = 1:NLupdates
                    [Lchain1(now),temp] = updateL_BD_3OUI1(E,LT1,lambda0(now),mu(now),sigma(now),lambda1(now),jumpPar,gamma_now,LT2,lambda2(now),negJ);                    
                    LT1 = Lchain1{now};
                    accRate(4) = accRate(4) + temp;
                    end
                case 2 % displacement
                    [Lchain1(now), temp] = updateL_dis_3OUI1(E,LT1,lambda0(now),mu(now),sigma(now),lambda1(now),jumpPar,LT2,lambda2(now),negJ,gamma_now);
                    accRate(5) = accRate(5) + temp;
            end
            % update all jump sizes with RW
            for ii = 1:NLupdates
            Lchain1(now) = updateL_sizeBlocked_3OU(E,Lchain1{now},lambda0(now),mu(now),sigma(now),lambda1(now),jumpPar,LT2,lambda2(now),negJ);
            accRate(6) = accRate(6) + temp;
            end
            %------------ update L2
            jumpPar = [2 beta2(now)];
            gamma3_now = struct('eta',gamma2(now),'t0',t0(now),'period',period,'theta1',theta1(now),'theta2',theta2(now));
            p = randi(2);
            negJflip = [negJ(2) negJ(1)];
            switch p
                case 1 % birth and death  
                    for ii = 1:NLupdates
                    Lchain2(now) = updateL_BD_3OUI1(E,LT2,lambda0(now),mu(now),sigma(now),lambda2(now),jumpPar,gamma3_now,LT1,lambda1(now),negJflip);                    
                    LT2 = Lchain2{now};
                    accRate(7) = accRate(7) + temp;
                    end
                case 2 % displacement
                    Lchain2(now) = updateL_dis_3OUI1(E,LT2,lambda0(now),mu(now),sigma(now),lambda2(now),jumpPar,LT1,lambda1(now),negJflip,gamma3_now);
                    accRate(8) = accRate(8) + temp;
            end
            for ii = 1:NLupdates
            Lchain2(now) = updateL_sizeBlocked_3OU(E,Lchain2{now},lambda0(now),mu(now),sigma(now),lambda2(now),jumpPar,Lchain1{now},lambda1(now),negJflip);
            accRate(9) = accRate(9) + temp;
            end
            
            %intermediate results
            if ~isempty(printR)
                if (and(i>printR(k),i<=printR(k)+n2))
                    j = j+1;                    
                    parameters(j,:) = [mu(now) lambda0(now) sigma(now) lambda1(now) theta2(now) lambda2(now) gamma2(now)/Tmax beta1(now) beta2(now) theta1(now) t0(now)];
                end
            end
            if j == n2
                j = 0; % reset counter and print intermediate results
                fprintf('%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n', ...
                    'mu', 'lambda0', 'sigma', 'lambda1','eta1','lambda2', 'eta2','beta1', 'beta2','delta1','theta1')
                fprintf('%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', mean(parameters(:,:)))
                fprintf('%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', std(parameters(:,:)))
%                 disp([mean(parameters(:,:));std(parameters(:,:))])
                k = min(k+1,length(printR)); % get next k
                %plotR(parameters2(:,4),1,[3 3 1],'lambda_2')
            end
            % save elements for output chain
            if 0 == mod(i,j3)
                j2 = j2 + 1;                
                par(j2,:) = [mu(now) lambda0(now) sigma(now) lambda1(now) theta2(now) lambda2(now) gamma2(now) beta1(now) beta2(now) theta1(now) t0(now)];
                Lchain_1(j2) = Lchain1(now);
                Lchain_2(j2) = Lchain2(now);
            end
            % swap before and now states
            now = 1 + mod(i,2);
            before = 1 + mod(i+1,2);
        end
                
    end % end gibbs3OUI1-Sigma

% Algorithm for 3-OU-I1 with new prior for Etai
    function gibbs3OUI1Eta
        for i = 2:n
            if 0==mod(i,10000)
                fprintf('MCMC step %d out of %d\n',i,n);
            end
            
            LT1 = Lchain1{before};
            LT2 = Lchain2{before};
            XL = E - negJ(1)*getY2(lambda1(before),LT1,Tmax+1) - negJ(2)*getY2(lambda2(before),LT2,Tmax+1);                        
            
            %------------ update mu
            mu(now) = update_mu(XL, lambda0(before), sigma(before), hPar.Amu, hPar.Bmu);
            
            %------------ update sigma
            sigma(now) = update_sigma(XL, lambda0(before), mu(now), hPar.Asigma, hPar.Bsigma);
            
            %------------ update lambda0
            [lambda0(now), temp]= update_lambda0(XL,lambda0(before),mu(now),sigma(now), sigmaPropLambda0);
            accRate(1) = accRate(1) + temp;
            
            %------------ update lambda1
            [lambda1(now), temp] = update_lambda1_3OU(E,LT1,lambda1(before),...
                lambda0(now),mu(now),sigma(now), sigmaPropLambda1, lambda2(before), LT2,negJ);
            accRate(2) = accRate(2) +temp;
                        
            %------------ update lambda2
            [lambda2(now), temp] = update_lambda2_3OU(E,LT2,lambda2(before),...
                lambda0(now),mu(now),sigma(now), sigmaPropLambda2, lambda1(now), LT1,negJ);
            accRate(3) = accRate(3) + temp;
                       
            %------------ update beta1
            beta1(now) = update_jumpBeta(LT1,hPar.Bbeta1);
            
            %------------ update beta2
            beta2(now) = update_jumpBeta(LT2,hPar.Bbeta2);
            
            %------------ update eta1
            % gamma1(now) = update_gamma(LT1,hPar.Beta1);
            
            %------------ update eta2
%             gamma2(now) = update_gamma(LT2,hPar.Beta2);
            % prior sensitivity for eta2
            x = gamma2(before)/Tmax;
            y = exp(normrnd(log(x),0.05));            
            NLT2 = length(LT2(1,LT2(1,:)>0));
            if NLT2>0
                L_L =  exp(-(y-x)*Tmax)*(y/x)^(NLT2) * y/x *(y>0)*(x>0);
            if rand <= L_L                
                gamma2(now) = y*Tmax;
                accRate(13) = accRate(13) + 1;
            else                
                gamma2(now) = x*Tmax;
            end
            end
            
            %------------ update intensity function parameters
            Tau = LT1(2,LT1(1,:)>0);
            [theta1(now), ~,t0(now), temp] = updateIntensityFunction(Tmax,Tau,t0(before),period,theta1(before),theta2(before));
            accRate([10 12]) = accRate([10 12]) + temp([1 3]);
            % prior sensitivity for eta1            
            Tau = Tau(2:end);
            x = theta2(before);
            y = exp(normrnd(log(x),0.05));    
            ratex = likelihoodPP(Tmax,0,Tau,t0(now),period,theta1(now),x);
            ratey = likelihoodPP(Tmax,0,Tau,t0(now),period,theta1(now),y);    
            L =  exp(ratey-ratex) * y/x * (y>0)*(x>0);
            if rand <= L
                theta2(now) = y;
                accRate(11) = accRate(11) + 1;
            else
                theta2(now) = x;
            end
            
            %------------ update L1
            % prepare jump parameters [distribution parameter1 parameter2]           
            jumpPar = [2 beta1(now)];    
            gamma_now = struct('eta',0,'t0',t0(now),'period',period,'theta1',theta1(now),'theta2',theta2(now));
            p = randi(2);
            switch p
                case 1 % birth and death  
                    for ii = 1:NLupdates
                    [Lchain1(now), temp] = updateL_BD_3OUI1(E,LT1,lambda0(now),mu(now),sigma(now),lambda1(now),jumpPar,gamma_now,LT2,lambda2(now),negJ);                    
                    LT1 = Lchain1{now};
                    accRate(4) = accRate(4) + temp;
                    end
                case 2 % displacement
                    [Lchain1(now),temp] = updateL_dis_3OUI1(E,LT1,lambda0(now),mu(now),sigma(now),lambda1(now),jumpPar,LT2,lambda2(now),negJ,gamma_now);
                    accRate(5) = accRate(5) + temp;
            end
            % update all jump sizes with RW
            for ii = 1:NLupdates
            [Lchain1(now),temp] = updateL_sizeBlocked_3OU(E,Lchain1{now},lambda0(now),mu(now),sigma(now),lambda1(now),jumpPar,LT2,lambda2(now),negJ);
            accRate(6) = accRate(6) + temp;
            end
            %------------ update L2
            jumpPar = [2 beta2(now)];
            gamma3_now = struct('eta',gamma2(now),'t0',t0(now),'period',period,'theta1',theta1(now),'theta2',theta2(now));
            p = randi(2);%p=0;
            negJflip = [negJ(2) negJ(1)];
            switch p
                case 1 % birth and death  
                    for ii = 1:NLupdates
                    [Lchain2(now), temp] = updateL_BD_3OUI1(E,LT2,lambda0(now),mu(now),sigma(now),lambda2(now),jumpPar,gamma3_now,LT1,lambda1(now),negJflip);                    
                    LT2 = Lchain2{now};
                    accRate(7) = accRate(7) + temp;
                    end
                case 2 % displacement
                    Lchain2(now) = updateL_dis_3OUI1(E,LT2,lambda0(now),mu(now),sigma(now),lambda2(now),jumpPar,LT1,lambda1(now),negJflip,gamma3_now);
                    accRate(8) = accRate(8) + temp;
            end
            for ii = 1:NLupdates
            Lchain2(now) = updateL_sizeBlocked_3OU(E,Lchain2{now},lambda0(now),mu(now),sigma(now),lambda2(now),jumpPar,Lchain1{now},lambda1(now),negJflip);
            accRate(9) = accRate(9) + temp;
            end
            
            %intermediate results
            if ~isempty(printR)
                if (and(i>printR(k),i<=printR(k)+n2))
                    j = j+1;                    
                    parameters(j,:) = [mu(now) lambda0(now) sigma(now) lambda1(now) theta2(now) lambda2(now) gamma2(now)/Tmax beta1(now) beta2(now) theta1(now) t0(now)];
                end
            end
            if j == n2
                j = 0; % reset counter and print intermediate results
                fprintf('%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n', ...
                    'mu', 'lambda0', 'sigma', 'lambda1','eta1','lambda2', 'eta2','beta1', 'beta2','delta1','theta1')
                fprintf('%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', mean(parameters(:,:)))
                fprintf('%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', std(parameters(:,:)))
%                 disp([mean(parameters(:,:));std(parameters(:,:))])
                k = min(k+1,length(printR)); % get next k
                %plotR(parameters2(:,4),1,[3 3 1],'lambda_2')
            end
            % save elements for output chain
            if 0 == mod(i,j3)
                j2 = j2 + 1;                
                par(j2,:) = [mu(now) lambda0(now) sigma(now) lambda1(now) theta2(now) lambda2(now) gamma2(now) beta1(now) beta2(now) theta1(now) t0(now)];
                Lchain_1(j2) = Lchain1(now);
                Lchain_2(j2) = Lchain2(now);
            end
            % swap before and now states
            now = 1 + mod(i,2);
            before = 1 + mod(i+1,2);
        end
                
    end % end gibbs3OUI1Eta

end % end gibbs function