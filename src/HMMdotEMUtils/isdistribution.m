function isP = isdistribution(candidate,rowOrCol)
% Determine if rows or columns of a matrix are probability distributions
% isP = isdistribution(candidate,rowOrCol)
% 
% rowOrCol is one of 'row'/1,'col'/2
% isP has same size as  sum(candidate,IND)
% where ind is 2 if row and 1 if col
%
% also see todistribution.m

if ischar(rowOrCol); rowOrCol=lower(rowOrCol); end
switch rowOrCol
  case {'row' 1}
    IND = 2;
  case {'col' 'column' 2}
    IND = 1;
  otherwise
    error('Second argument should be ''row''/1 or ''col''/2');
end
sumTol = size(candidate,IND)*eps;

isP = abs(sum(candidate,IND)-1)<sumTol;
end
