classdef FunMinimizer < handle
% Helps pass arguments to functions that you might use 
% to update the conditional parameters of the model
%
% This class is a handle subclass and so acts like a pointer

properties(SetAccess=private)
  % The minimizer's name
  % A valid minimizer takes in arguments:
  %   "funToMinimize" - a function that is to be minimized
  %   "x0" - something to be minimized over
  %   "options" - some form of options for the function
  % 
  % e.g. should act like 'fminsearch'
  funName = '';
end
properties
   % Params that can be passed to funName
   % e.g. optimset('fminsearch')
  funParams = struct;
end

methods
  function obj = FunMinimizer(funName)
  % Construct an FunMinimizer by specifying a string
  % obj = FunMinimizer(funName)
  % Supported are: 'fminsearch'
  if nargin>0
    switch funName
      case 'fminsearch'
        obj.funName = funName;
        obj.funParams = optimset(obj.funName);
      otherwise
        error(['function ' funName ' not supported by FunMinimizer']);
    end
  end
  end
  
  function result = feval(obj,funToMinimize,x0)
    % result = feval(obj,funToMinimize,x0)
    %
    % feval for this object
    result = feval(obj.funName,funToMinimize,x0,obj.funParams);
  end
end
end
