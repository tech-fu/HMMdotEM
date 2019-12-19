function [P S] = todistribution(A,rowOrCol)

% Divides each column of a matrix by its sum.
% [P S] = tocoldistribution(A,rowOrCol)
%
% rowOrCol can be 'row'/1 ,'col'/2/'column'
% Optionally returns the sums S

if ischar(rowOrCol); rowOrCol=lower(rowOrCol); end
switch lower(rowOrCol)
  case {'row' 1}
    S = sum(A,2);
  case {'col' 'column' 2}
    S = sum(A,1);
  otherwise
    error('Second argument should be ''row''/1 or ''col''/2');
end
if issparse(A)
  L = length(S); % Cannot be 3D
  P = sparse(1:L,1:L,1./S)*A;
else
  P = bsxfun(@ldivide,S,A);
end
