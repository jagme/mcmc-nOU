function y = update_gamma(L, Bgamma)
%UPDATE_GAMMA Gibbs sampler for gamma = eta*Tmax,
% eta is the jump rate of the compound Poisson process L
%
% Input arguments:
%
%   L - jump process L as 2-by-p matrix. The first row contains the jump
%   sizes, whilst the the second row contains the time grid
%
%   Bgamma - prior hyperparameter
%
% Output arguments:
%
%   y - sample of gamma

Tmax = max(L(2,:));
idx = find(L(1,:)>0);
sum_jumps = length(idx);

A = 1;
B = Bgamma;

y = gamrnd(A+sum_jumps , 1/(B+Tmax),1,1);
y = y*Tmax;
