function [Y22, Y2] = getY2(lambdai, L, nOut)
%GETY2 Compute Non-Gaussian OU process Yi, i>0 given the background 
%compound Poisson process L, lambdai and initial condition Yi(0) = 0.
% The process is returned at times 
% t0 = 0, t1 = 1, t2 = 2, ..., tn = T
%
% [Y22, Y2] = getY2(lamb_2,L,n_out)
%
% Input arguments: 
%
%   lambdai - inverse of speed of mean reversion of process Yi
%
%   L - 2-by-p matrix representing the compound Poisson process Li of Yi.
%       The first row contains the jump sizes at the times given in the second
%       row. p is variable and depends on the total number of jumps on L
%
%   nOUt - number of observations of the data X, ie n + 1
%
% Output arguments:
%
%   Y22 - the process Yi at times 0, 1, ..., T. Y22 is of length nOut
%
%   Y2 - the exact path of Yi implied by L, lambdai and initial condition
%   Yi(0) = 0
%

% To improve computational performance, this routine has been compiled 
% into a C++ Mex file and is available for Unix and Linux platforms with 
% Matlab 2014b or later. Windows platforms will use this 
% pure Matlab implementation automatically. If calling this function under
% Unix/Linux produces an error, probably your platform/version is not 
% compatible.  A workaround is to rename the files
% getY2.mexmaci64 (Unix) and/or getY2.mexa64 (Linux), e.g.
% getY2.mexmaci64 ->>> getY2_.mexmaci64 so the mex files are not called

% %{
lenT = length(L(2,:));  % length of process L

Y2  = zeros(1,lenT);    % exact OU process based on L
Y22 = zeros(1, nOut);   % over fix dates

dt1  = L(2,2:end) - L(2,1:end-1);
temp = exp(-1/lambdai*dt1);
k = 1;
for j = 1:lenT-1    
    Y2(j+1) = temp(j)*Y2(j) + L(1,j+1);
    if L(2,j+1) == k
        Y22(k+1) = Y2(j+1); % get only fix dates 1, 2, ...
        k = k + 1;
    end
end

%}
