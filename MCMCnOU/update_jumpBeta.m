function y = update_jumpBeta(L, B)
%UPDATE_JUMPBETA Gibss sampler for mean jump size beta of Exp. distribution
%
% Input arguments:
%
%   L - jump process L as 2-by-p matrix. The first row contains the jump
%   sizes, whilst the the second row contains the time grid.
%
% Output argument:
%
% y - sample of beta
%

idx = find(L(1,:)>0);   % get index of Exponential observations (jump sizes)
n = length(idx);        % number of observations

%prior parameters for 1/beta, Gamma(A,B) mean A/B
A = 1;
% B = 1;

y = gamrnd(A+n,1/(B+sum(L(1,idx))),1,1);
y = 1/y; % return Inverse-Gamma sample
