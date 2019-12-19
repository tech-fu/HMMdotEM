function x = randgauss(m,S,numCols)
% Generate random gaussian with mean and standard deviation matrix
% x = randgauss(m,S,numCols)
% The point is to avoid the Statistics toolbox dependency
%
% varargin is anything that can be passed to randn
% m - mean column vector
% S - standard deviation matrix
%
% x = bsxfun(@plus,S*randn(length(m),numCols),m);

assert(hasonlydim(m,1),'mean must be a column vector');
x = bsxfun(@plus,S*randn(length(m),numCols),m);

end
