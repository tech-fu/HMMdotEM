function bool = hasonlydim(v,dim)
% Determine if a tensor has size>1 only for dimension dim
% bool = hasonlydim(v,dim)
%
% bool = size(v,dim) == numel(v);
%
% Also see iscol

bool = size(v,dim) == numel(v);

