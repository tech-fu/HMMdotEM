classdef HMMParams < handle
% A class that stores HMM params and inherits from handle so sharing is easier.
% Has properties that HMM model has (only the properties that can change)
% Can be given to multiple models and when one model is changed
% all other models will also change.
% This allows easy handling of data that comes in segments with
% almost no modification to existing code.
%
% No constructor is given, just call HMMParams() to instantiate

properties(SetObservable) % SetObservable allows event on change
  startP % param that should be sharable
  transitionP % param that should be sharable
  condParams % param that should be sharable
  condParamNames % param that should be sharable
  % numHids and condLogLikFun are not necessary since they are private methods
end

end
