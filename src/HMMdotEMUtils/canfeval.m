function canEval = canfeval(candidate)
% Determing whether a variable can be used in feval
%
% canEval = canfeval(candidate)
% candidate - any variable
% canEval - true iff candidate is a row string or function_handle
%           or has a method called "feval"

canEval = isa(candidate,'function_handle') ||...
  (ischar(candidate) && hasonlydim(candidate,2)) ||... % Can be string
  any(strcmp(methods(candidate),'feval'));
