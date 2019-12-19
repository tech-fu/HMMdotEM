classdef HMMTemplate < HMM
% A template that can be reused to define subclasses of HMM for specific models
% Includes commonly redefined methods
% IMPORTANT: you MUST DEFINE the methods in the first method block,
%   such as "condloglikfun" and the constructor.
% Replace all occurences of "HMMTemplate" with whatever you want to call your
% model/subclass (e.g. "MyHMM").
%
% See the original HMM class if methods are not self-explanatory
% If you want to use the original HMM methods don't include/rewrite them in you
% subclass, only include the methods you would like to redefine,
% this means you should completely remove any code that you don't want to use
% (except of course the code that is mandatory).

%% THESE MUST BE IMPLEMENTED
methods 
  function obj = HMMTemplate(varargin)
    % This one needs to be redefined
    % Change "HMMTemplate" above to whatever your class is called
    obj = obj@HMM(varargin{:});
  end
  
  function condLoglLik = condloglikfun(obj,condParamsArg,dataArg)
    % This is the meat of the class
    % Checking the inputs is not necessary, but may be helpful
    assert(all(size(condParamsArg)==[obj.numCondParams obj.numHids]),...
      'condParamsArg size mismatch');
    assert(size(dataArg,1)==obj.dataDim,...
      'dataArg size did not match expected');
  end
end

%% These are optional.
% Completely remove them is you don't need to alter them
methods
  function obj = initcondparams(obj,numCondParams)
    % Initialize the conditional parameters given just how many to use
  end
  
  function [newParams oldVal newVal] = fitcondparams(obj,FM)
    % This is the function that is used in the mstep and
    % can be used to provide an analytic solution to the data term minimization
    % The first arg will be a FunMinimizer object
    % oldVal and newVal are the data term values before and after the fit
    % their difference signifies the gain from the data change in condParams
    assert(obj.isparamgroup,'This HMM is not a paramgroup');
    respAll = [obj.computeresp.resp];
    dataAll = [obj.data];
    oldVal = 0;
    newVal = 0;
  end
  
  function obj = namecondparams(obj)
    % Names for conditional parameters are stored in obj.Params.condParamNames
    obj.test(obj.isparamgroup,'obj is not a param group');
  end
  
  function [data hidPath] = sampledata(obj,numData)
    dDim = obj.dataDim;
    numH = obj.numHids;
    hidPath = obj.samplepath(numData);
    data = zeros(dDim,numData);
    isHid = accumarray(hidPath,1:numData,[numH 1],@(x) {x},{[]});
    numOccur = cellfun(@length,isHid);
    for k=1:obj.numH
      % Make all the samples ones
      % MODIFY THE LINES BELOW (NOW SAMPLES ALL nans)
      data(:,isHid{k}) = nan(dDim,numOccur(k));
    end
  end
  
  function displaydata(obj,titlePrefix,useViterbi)
    assert(isscalar(obj),'cannot display data for an array of HMMs');
    if nargin<2
      titlePrefix = '';
    end
    if nargin<3
      useViterbi = true;
    end
    
    H = obj.numHids;
    allMarkers = allmarkers;
    numMarkers = size(allMarkers,1);
    needMarkers = H-numMarkers;
    if needMarkers>0 % If the markers are not enough pad with one value
      allMarkers(numMarkers:numMarkers+needMarkers,:) = 'k.';
    end
    markSize = obj.displaymarkersize;
    
    % Insert the rest of the code here
    
  end
end

end
