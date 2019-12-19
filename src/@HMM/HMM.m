classdef HMM
% HMM model + parameters are in this class
% The main properties are the data, the number of hidden variables
% and the arbitrary log-likelihood function.
% See the property descriptions for details
%
% IMPORTANT: This class has (some) Abstract methods
%   (search for "methods(Abstract)" to find them).
%   To be able to use the class you must define the subclass 
%   that implements your model. Don't worry most of the work is done already...
%   See HMMTemplate.m and GaussMixHMM.m for examples.
%
% Some notation helpful in comments
% X_t,Y_t - data t \in 1:T (Starts at time-step #1)
% Z_t   - hidden variables multinomial 1:K
% a_k = P(Z_1=k)
% \theta_k - P(X_t,Y_t | Z_t=k) ~ F_{\theta_k}
% A_{ij} = P(Z_{t+1} = j | Z_t = i)

%% Constants
properties(Constant)
  % Quantity to add to transition probs to guarantee process is ergodic
  % Will only be used if the property assumeEquilibrium is true
  ERGODICMIN = .01;
end

%% Immutable properties (Changing these will mean a different model)
properties(SetAccess=protected)
  data % dataDim x numTimeSteps data matrix
  numHids % Number of hidden variables (K)
  
end

%% Important properties that should be changed only by class methods
properties(SetAccess=private,GetAccess=private)
  % Contains the parameters and facilitates sharing among HMM's
  % A key component in the training of a model with more than one
  % set of data.
  % See the dependent properties of this class.
  % See class HMMParams
  Params 
  
end

%% Parameters that can be set as you like
properties
  %% Assumptions
  
  % Assume the data is coming from the equilibrium distribution?
  % Will effectively use the left 1-eigenvector of transitionP as startP
  % So when you call obj.startP it will return the eigenvector
  % Note: regardless of this obj.setstartp(startP) will set the property startP
  % even though it is not returned by obj.startP
  assumeEquilibrium = true;
  
  %% Flags
  
  verbose = false; % Print stuff?
  display = false; % Display stuff? (e.g. Figures)
  testing = false; % Test code using assertions, switch off to speed up code
end

%% Dependent quantities
properties(Dependent)
  %%
  % Dimension of the data ( (X,Y) has dimension 2 )
  dataDim;
  
  %% Parameters from Params property
  
  condParams % numCondParams x numHids params (condParams(:,k) is \theta_k)
  condParamNames % You can assign names to the parameters (for display)
  
  startP % 1 x numHids probabilities (startP(k) is a_k)
  
  % numHids x numHids transition probabilities (rows sum to 1)
  % (transitionP(i,j) is A_{ij})
  transitionP; 
  
  % numFreeParams is numHids*(numHids-1) + numHids*numCondParams
  %           + (numHids-1) (start probabilities are added if not equilibrium)
  numFreeParams;
  
  %% Sizes
  
  numTimeSteps % The time dimension of the data
  
  numCondParams % The number of parameters each state has (the size of \theta_k)
end

%% Temporary quantities dependent on model parameters
% The method emptyparamdependent empties these
properties(Transient,SetAccess=protected)
  %% Forward-back
  
  condLogLik % Obtained directly from condloglikfun log(P(data_t | Z_t))
  condLik % P(data_t | Z_t). Computed and stored to avoid re-exponentiation
  
  alphaS % P(Z_t | data_1:t) 
  scaling  % P(data_t | data_1:{t-1})
  hmmLogLik % P(data) - data log-likelihood
  
  betaS  % P(data_{t+1}:T | Z_t)/P(data_{t+1}:T | data_1:t)
  
  resp % Responsibilities for each time step P(Z_t | data)
  transitionStats % \sum_{t>1} P(Z_{t+1} , Z_t | data)
  
  viterbiPath % Integer vector containing the most probable path {Z_t}_max
  viterbiLogLik % P({Z_t}_max,data)
  % IMPORTANT!!! If adding another property also updata "emptyparamdependent"
end

methods(Abstract) % Each HMM subclass has to define these.
  condloglikfun(obj,condParams,data)
end

%% CONSTRUCTOR and INITIALIZATION
methods
  function obj = HMM(data,numHids,condParamArg)
    % obj = HMM(data,numHids,condParams)
    % data is DATADIMxNUMTIMESTEPS
    % numHids - scalar, natural
    % condParams - initial parameters of conditional distribution
    %          
    % obj = HMM(data,numHids,numCondParams)
    % numCondParams - number of conditional parameters per state
    %                 The dimension of \theta_k
    %   The condParams will be estimated using the method 'initcondparams'
    
    if nargin>0 % A matlab Object-Oriented issue
      %% Assign properties
      obj = obj.setdata(data);
      
      assert(isscalar(numHids),'numHids should be a scalar value');
      obj.numHids = numHids;
      obj.Params = HMMParams();
      obj = obj.initializeP;
      if isscalar(condParamArg) % Must be 
        numCondParams = condParamArg;
        obj = obj.initcondparams(numCondParams);
      else % Must be matrix of parameters
        assert(~isempty(condParamArg) && size(condParamArg,2)==obj.numHids,...
        'condParams should be numCondParams x numHids matrix');
        obj.Params.condParams = condParamArg;
      end
      
      obj.namecondparams();
      
      %% Check if the function can be evaluated
      obj.checkcondloglikfun; % Does it return the right form?
    end
  end
  function objArray = constructparamgroup(obj,dataCell)
    % Construct a "param group" of HMM's that have a common set of parameters
    % objArray = obj.constructhmmgroup(dataCell)
    % dataCell is a cell column containing data in each cell
    %
    % This is the same as the regular constructor but it returns and array
    % of HMM's which have the same Params handle instance
    % The group has all properties equal except for data.
    % Modifying obj(1).transitionP will also change obj(2).transitionP
    % so these are gauranteed to have the same parameters forever.
    %
    % the assumeEquilibrium flag will be set to true for each HMM in the array
    % you can change this for each one individually after (it is not shared)
    assert(hasonlydim(dataCell,1),'dataCell should be column');
    numSegments = length(dataCell);
    objArray(numSegments,1) = obj;
    for h = 1:numSegments
      objArray(h,1) = objArray(numSegments,1).setdata(dataCell{h});
    end
    [objArray.assumeEquilibrium] = deal(true);
    assert(objArray.isparamgroup,'Failed to construct param group of HMM''s');
  end
  
  %% Initialize
  function obj = initializeP(obj)
    % Sets the start and transition probabilities to uniform
    % obj = obj.initializeP()
    obj = obj.setstartp(todistribution(ones(1,obj.numHids),'row'));
    obj = obj.settransitionp(todistribution(ones(obj.numHids),'row'));
  end
  
  function obj = initcondparams(obj,numCondParams)
    % Estimate/initialize conditional parameters (not based on model)
    % obj = obj.initcondparams(numCondParams)
    % 
    % Initialization is very hard for arbitrary params
    % but this methods can be redefined in subclasses
    % This one just uses rand(numCondParams,obj.numHids)
    warning('HMM:initcondparams:poor',...
      'This general version of initcondparams method is very poor');
    obj = obj.setcondparams(rand(numCondParams,obj.numHids));
  end
  
  %% Check
  function obj = checkcondloglikfun(obj)
    % Check if the conditional likelihood function returns proper result
    % when evaluated on the whole data set for each state
    % obj.checkcondloglikfun
    %
    % Note: canfeval.m only checks if feval can be called, but not the output
    % this checks only the output.
    testingStore = obj.testing;
    obj.testing = false; % Avoid checks in computecondloglik
    obj = obj.computecondloglik;
    logL = obj.condLogLik;
    assert(all(size(logL)==[obj.numHids obj.numTimeSteps]) ...
        && isreal(logL),...
      'The conditional log-likelihood function didn''t return the proper result')
    obj.testing = testingStore;
  end
end

%% SET/GET
methods
  %% Set
  function obj = setdata(obj,data)
    % Changes data property (use with extreme caution)
    % obj = obj.setdata(dataMat)
    numDims = obj.dataDim;
    assert(numDims==0 || size(data,1)==numDims,...
        ['data should be ' num2str(numDims) 'x numTimeSteps']);
    obj.data = data;
    obj = obj.emptyparamdependent;
  end
  
  function obj = setstartp(obj,sP)
    % Set the start probability vector
    % obj = obj.setstartp(sP)
    % Can deal with param groups of HMMs
    obj.testparamgroup;
    assert(all(size(sP)==[1 obj(1).numHids]),...
      'startP can only be set to 1xnumHids')
    assert(isdistribution(sP,'row'),...
      'cannot set startP to non-normalized distribution (see todistribution.m)')
    obj(1).Params.startP = sP;
    obj = obj.emptyparamdependent;
  end
  
  function obj = settransitionp(obj,tP)
    % Set the transition probability matrix
    % obj = obj.settransitionp(tP)
    % Can deal with param groups of HMMs
    % The distribution will be fixed if process is not ergodic and 
    % the property assumeEquilibrium is true
    obj.testparamgroup;
    assert(all(size(tP)==obj(1).numHids*[1 1]),...
      'transitionP can only be set to numHidsxnumHids')
    assert(all(isdistribution(tP,'row')),...
      'cannot set transitionP to non-normalized distribution (see todistribution.m)')
    if any([obj.assumeEquilibrium])
      tP = todistribution(max(tP,HMM.ERGODICMIN),'row');
    end
    obj(1).Params.transitionP = tP;
    obj = obj.emptyparamdependent;
  end
  
  function obj = setcondparams(obj,condParams)
    % Set the conditional parameter matrix
    % obj = obj.setcondparams(condParams)
    % Can deal with param groups of HMMs
    obj.testparamgroup;
    assert(size(condParams,2)==obj(1).numHids,...
      'condParams is not ? x numHids');
    assert(isempty(obj.condParams) ||...
      all(size(condParams)==size(obj.condParams)),...
      'condParams is trying to change size');
    obj(1).Params.condParams = condParams;
    obj = obj.emptyparamdependent;
  end
  function obj = namecondparams(obj,nameCell)
    obj.testparamgroup;
    if nargin<2
      nameCell = cellfun(@num2str,num2cell(1:obj.numCondParams),...
                                                'UniformOutput',false);
    else
      assert(length(nameCell)==obj(1).numCondParams,...
        'cond param names (nameCell) are not right length');
      isrowchar = @(s) hasonlydim(s,2) && ischar(s);
      assert(all(cellfun(isrowchar,nameCell)),...
        'nameCell should contain row char');
    end
    obj(1).Params.condParamNames = nameCell;
  end
  
  %% Get
  function dataDim = get.dataDim(obj)
    dataDim = size(obj.data,1);
  end
  function condParams = get.condParams(obj)
    condParams = obj.Params.condParams;
  end
  function nameCell = get.condParamNames(obj)
    nameCell = obj.Params.condParamNames;
  end
  function startP = get.startP(obj)
    % Returns the starting probability
    % If assumeEquilibrium will return equilibrium distribution
    if obj.assumeEquilibrium
      startP = obj.computeequilibrium;
    else
      startP = obj.Params.startP;
    end
  end
  function transitionP = get.transitionP(obj)
    transitionP = obj.Params.transitionP;
  end
  
  function numTimeSteps = get.numTimeSteps(obj)
    numTimeSteps = size(obj.data,2);
  end
  
  function numCondParams = get.numCondParams(obj)
    numCondParams = size(obj.condParams,1);
  end
  
  function nFP = get.numFreeParams(obj)
    % numFreeParams is numHids*(numHids-1) + numHids*numCondParams
    %           + numHids (start probabilities are added if not equilibrium)
    nFP = (1-obj.assumeEquilibrium + obj.numHids)*(obj.numHids-1) + ...
                                                obj.numHids*obj.numCondParams;
  end
  
  function bool = isparamgroup(objArray)
    % Check if as array of HMM's is of type "param group"
    % bool = objArray.isparamgroup;
    % Basically they must have the same main properties except for data.
    % Most importantly the Params property should be the same handle
    % To be a param group the array must be column
    % Scalar objArrays are trivially param groups
    bool = hasonlydim(objArray,1);
    i=2;
    while bool && i<=length(objArray)
      bool = objArray(1).Params == objArray(i).Params &&... 
             objArray(1).numHids == objArray(i).numHids;
      i = i+1;
    end
  end
  function testparamgroup(obj)
    if any([obj.testing])
      assert(obj.isparamgroup,'obj is not a param group. Need to loop.');
    end
  end
  
  function logL = paramgrouploglik(objArray)
    % Compute a column of log-likelihoods for an HMM param group
    % logL = objArray.fullparamgrouploglik
    %
    % logL is length(objArray)x1
    numObj = length(objArray);
    logL = -inf(numObj,1);
    for i=1:numObj
      obj = objArray(i).computehmmloglik;
      logL(i) = obj.hmmLogLik;
    end
  end
  
end

%% Compute Model Quantities
methods
  %% Conditional log-L
  function objArray = computecondloglik(objArray)
    % Update the properties "condLogLik"
    % obj = obj.computecondloglik()
    % obj.condLogLik is obtained from the method condloglikfun
    
    for i=1:numel(objArray)
      obj = objArray(i);
      if isempty(obj.condLogLik)
        obj.condLogLik = obj.condloglikfun(obj.condParams,obj.data);
        if obj.testing
          assert(~any(isnan(obj.condLogLik(:))),'Found condLogLik nans');
          assert( isreal(obj.condLogLik),'Found condLogLik complex');
        end
      end
      objArray(i) = obj;
    end
  end
  function objArray = computecondlik(objArray)
    % obj.condLik contains the exp of "condLogLik"
    for i=1:numel(objArray)
      obj = objArray(i);
      if isempty(obj.condLik)
        obj = obj.computecondloglik;
        obj.condLik = exp(obj.condLogLik);
        if obj.testing
          assert(~any(all(obj.condLik==0,1)),...
            'Found all 0 columns in condLik. P(data)=0. Must fix parameters.');
        end
      end
      objArray(i) = obj;
    end
  end
  
  %% Forward-back
  function objArray = forwardpass(objArray)
    % Update the quantities alphaS,scaling
    % obj = obj.forwardpass()
    for i=1:numel(objArray)
      obj = objArray(i);
      sP = obj.startP;
      tP = obj.transitionP;
      if isempty(obj.alphaS) || isempty(obj.scaling)
        obj = obj.computecondlik; % Should not recompute existing quantities
        [obj.alphaS obj.scaling] = HMM.alphahelper(obj.condLik,sP,tP);
        if obj.testing
          assert(~any(isnan(obj.alphaS(:))),'Found alphaS nans');
          assert(~any(isnan(obj.scaling)),'Found scaling nans');
          assert(~any(obj.scaling==0),'Found 0 scaling');
        end
      end
      objArray(i) = obj;
    end
  end
  function objArray = computehmmloglik(objArray)
    % Update the property hmmLogLik
    % obj = obj.computehmmloglik()
    for i=1:numel(objArray)
      obj = objArray(i);
      if isempty(obj.hmmLogLik)
        obj = obj.forwardpass; % Should not recompute existing quantities
        obj.hmmLogLik = sum(log(obj.scaling),2);
        if obj.testing
          assert(isscalar(obj.hmmLogLik),'hmmLogLik is not a scalar');
          assert(~isnan(obj.hmmLogLik) && isreal(obj.hmmLogLik),...
          'hmmLogLik nan or complex');
        end
      end
      objArray(i) = obj;
    end
  end
  function objArray = backpass(objArray)
    % Update the property "betaS"
    % obj = obj.backpass()
    for i=1:numel(objArray)
      obj = objArray(i);
      if isempty(obj.betaS) % Don't recompute if already have it
        obj = obj.forwardpass; % Should not recompute existing quantities
        tP = obj.transitionP;
        obj.betaS = HMM.betahelper(obj.condLik,obj.alphaS,obj.scaling,tP);
        if obj.testing
          assert(~any(isnan(obj.betaS(:))),'Found betaS nans');
        end
      end
      objArray(i) = obj;
    end
  end
  function objArray = forwardback(objArray)
    % Complete the forward and back passes
    % obj = obj.forwardback()
    for i=1:numel(objArray)
      obj = objArray(i);
      obj = obj.forwardpass; % This is not necessary, but doesn't matter
      obj = obj.backpass;
      objArray(i) = obj;
    end
  end
  
  %% Other quantities
  function objArray = computeresp(objArray)
    % Update the property "resp"
    for i=1:numel(objArray)
      obj = objArray(i);
      if isempty(obj.resp)
        obj = obj.forwardback; % Should not recompute existing quantities
        obj.resp = obj.alphaS.*obj.betaS;
        if obj.testing
          assert(~any(all(obj.resp==0,1)),...
            'Found t such that P(Z_t=k|data) = 0 for all k.');
        end
      end
      objArray(i) = obj;
    end
  end
  function objArray = computetransitionstats(objArray)
    % Update the property "transitionStats"
    % obj = obj.computetransitionstats();
    for i=1:numel(objArray)
      obj = objArray(i);
      if isempty(obj.transitionStats)
        obj = obj.forwardback; % Should not recompute existing quantities
        tP = obj.transitionP;
        temp = bsxfun(@rdivide, obj.condLik(:,2:end) , obj.scaling(1,2:end));
        temp = temp .* obj.betaS(:,2:end);
        obj.transitionStats = (obj.alphaS(:,1:end-1)*temp').*tP;
      end
      objArray(i) = obj;
    end
  end
  
  function objArray = computeviterbi(objArray)
    % Update the properties viterbiPath and viterbiLogLik
    % obj = obj.computeviterbi()
    for i=1:numel(objArray)
      obj = objArray(i);
      if isempty(obj.viterbiPath) || isempty(obj.viterbiLogLik)
        obj = obj.computecondloglik;
        logTransP = log(obj.transitionP);
        logStartP = log(obj.startP);
        [obj.viterbiPath obj.viterbiLogLik] =...
          HMM.viterbihelper(obj.condLogLik,logStartP,logTransP);
      end
      objArray(i) = obj;
    end
  end
  
  %% EM-usefuls
  function [newParams oldVal newVal] = fitcondparams(obj,FM)
    % [newParams dataTermDelta] = obj.fitcondparams(FM)
    % FM - FunMinimizer object or similar
    %      has to be a MINIMIZER!
    % 
    % newParams - are the condParams that MAXIMIZE the data term in the M-step
    % delta - is the increase in the data term from the result of this function
    obj.testparamgroup;
    respAll = [obj.resp]; % Works only for row-time data
    dataAll = [obj.data];
    
    qFun = @(cp) -sum(sum(respAll.*obj(1).condloglikfun(cp,dataAll)));
    oldVal = qFun(obj(1).condParams);
    newParams = feval(FM,qFun,obj(1).condParams);
    newVal = qFun(newParams);
  end
end

%% Methods for Memory Usage (VERY IMPORTANT!)
methods
  function obj = emptyforwardback(obj)
    % Empty condLogLik,alphaS,scaling,betaS
    % Can be used in EM to reduce memory
    % Can deal with arrays of HMMs
    [obj.condLogLik] = deal([]);
    [obj.condLik] = deal([]);
    [obj.alphaS] = deal([]);
    [obj.scaling] = deal([]);
    [obj.betaS] = deal([]);
  end
  function obj = emptyparamdependent(obj)
    % Empty all the properties that depend on model parameters
    % obj.emptyparamdependent
    % Can deal with arrays of HMMs
    obj = obj.emptyforwardback;
    [obj.hmmLogLik] = deal([]);
    [obj.resp] = deal([]);
    [obj.transitionStats] = deal([]);
    [obj.viterbiPath] = deal([]);
    [obj.viterbiLogLik] = deal([]);
  end
end

%% Sampling/MC
methods
  function hidPath = samplepath(obj,numTimeSteps)
    % Sample a path of hidden states from the model
    % hidPath = obj.samplepath(numTimeSteps)
    % hidPath is 1 x numTimeSteps
    %  hidPath(t)==k is true iff state k was sample at time t
    
    hidPath = obj.samplepathhelper(numTimeSteps,obj.startP,obj.transitionP);
  end
  
  function p = computeequilibrium(obj)
    % Compute the equilibrium distribution over states
    % p = obj.computeequilibrium
    [V D] = eig(obj.transitionP'); % Not necessarily sorted.
    [d dI] = min(abs(diag(D) - 1));
    p = todistribution(V(:,dI(1))','row');
    if obj.testing
      assert(obj.numHids==1 || all(p>0 & p<1),...
        'Found <=0 or >=1 in equilibrium distribution');
    end
  end
end

%% Model Selection
methods
  function logLik = testloglik(obj,testData)
    TestModel = obj.setdata(testData).computehmmloglik;
    logLik = TestModel.hmmLogLik;
  end
  function s = bic(obj)
    % Compute the BIC (Bayesian Information Criterion) for the model
    % s = obj.bic
    % s is hmmloglik - .5 * numFreeParams * log(numTimeSteps)
    obj = obj.computehmmloglik;
    s = obj.hmmLogLik - (obj.numFreeParams*log(obj.numTimeSteps))/2;
  end
  function s = icl(obj)
    % Compute the ICL model selection criterion for the model
    % s = obj.icl
    % s is based on the likelihood of the most probable path
    obj = obj.computeviterbi;
    s = obj.viterbiLogLik - (obj.numFreeParams*log(obj.numTimeSteps))/2;
  end
end

%% IO AND DISPLAY
methods
  function fprintf(obj,varargin)
    % Calls fprintf if obj.verbose
    % obj.fprintf(varargin)
    if obj.verbose
      fprintf(varargin{:});
    end
  end
  
  function markerSize = displaymarkersize(obj)
    T = obj.numTimeSteps;
    groupThres = 10.^(2:5);
    markerSize = 1+nnz(T<=groupThres);
  end
  
  function displayhmm(obj)
    % Plot startP, transitionP and condParams
    % obj.displayhmm()
    %
    % Uses subplots(1,3,:) respectively
    H = obj.numHids;
    
    subplot(1,3,1);hold on;
    sP = obj.startP(:);
    imagesc(sP,[0 1]);
    axis ij;
    for i=1:length(sP)
      text(1,i,num2str(sP(i),'%.2f'));
    end
    if obj.assumeEquilibrium
      s = 'Equilibrium distribution';
    else
      s = 'Start Probs';
    end
    title(s);
    xlabel('white=1,black=0'); ylabel('state')
    set(gca,'XLim',[.5 1.5],'YLim',[.5 H+.5]);
    
    subplot(1,3,2);hold on;
    tP = obj.transitionP;
    imagesc(tP,[0 1]); colormap(gray);
    axis ij;
    for i=1:H
      for j=1:H
      text(j,i,num2str(tP(i,j),'%.2f'));
      end
    end
    title('Transition Probs');
    xlabel('to state'); ylabel('from state');
    set(gca,'XLim',[.5 H+.5],'YLim',[.5 H+.5]);
    
    subplot(1,3,3);hold on
    obj.displaycondparams;
    
  end
  function displaycondparams(obj)
    % Plot the "condParams" matrix
    % obj.displaycondparams()
    H = obj.numHids;
    C = obj.numCondParams;
    for i=1:C
      for h=1:H
        text(h,i,num2str(obj.condParams(i,h),' %.2g'));
      end
    end
    set(gca,'XLim',[.5 H+.5],'YLim',[.5 C+.5],...
      'YTick',1:C,'XTick',1:H,'YTickLabel',obj.condParamNames);
    axis ij
    title('Cond. Params');
    xlabel('Hidden State'); ylabel('Parameter');
  end
  
  displaydata(obj,titlePrefix,useViterbi);
  
  displayhids(obj,titlePrefix,useViterbi);
end

%% STATIC
methods(Static)
  
  %% Helpers
  function [alphaS scaling] = alphahelper(condLik,sP,tP)
    % Compute the scaled alpha quantities and scaling
    % for a set of conditional likelihoods, start probs and transition probs
    % [alphaS scaling] = alphahelper(condLik,sP,tP)
    if exist('fastalphascaled','file')==3 % Use the mex file (mex/dll==3)
      [alphaS scaling] = fastalphascaled(condLik,sP,tP);
    else
    % Should not be too slow since the loop has only one statement in it
    % Try "feature accel on" for speed up.
      numTimeSteps = size(condLik,2);
      alphaS = nan(size(condLik));
      scaling = nan(1,numTimeSteps);
      [alphaS(:,1) scaling(:,1)] =...
                    todistribution(condLik(:,1).*sP(:),'col');
      % Recurse for alpha
      for t=2:numTimeSteps
        % P(Z_t | data_1:t) \propto 
        % P(data_t | Z_t) \sum_{k} P(Z_t|Z_{t-1}=k) P(Z_{t-1} | data_{1:t-1})
        [alphaS(:,t) scaling(:,t)] = todistribution(...
                                    condLik(:,t).*(tP'*alphaS(:,t-1)),...
                                                    'col');
      end
    end
  end
  function betaS = betahelper(condLik,alphaS,scaling,tP)
    % Compute the scaled beta quantities
    % for a set of conditional likelihoods, scaled alpha, scaling
    % and transition probs
    % betaS = betahelper(condLik,alphaS,scaling,tP)
    if exist('fastbetascaled','file')==3 % Use the mex file (mex/dll==3)
      betaS = fastbetascaled(condLik,alphaS,scaling,tP);
    else
      % Should not be too slow since the loop has only one statement in it
      % Try "feature accel on" for speed up.
      numTimeSteps = size(condLik,2);
      betaS = nan(size(condLik));
    
      betaS(:,end) = 1;
      for t=(numTimeSteps-1:-1:1)
        % P(data_{t+1}:T | Z_t)/P(data_{t+1}:T | data_1:t) =
        %  \sum_k P(data_{t+2}:T | Z_{t+1}=k) /
        %         P(data_{t+2}:T | data_1:{t+1}) 
        %    *    P(data_{t+1}|Z_{t+1}=k) P(Z_{t+1}=k|Z_t) / 
        %         P(data_{t+1ss}|data_1:t)

        betaS(:,t) = tP*( betaS(:,t+1).*( condLik(:,t+1)./scaling(1,t+1) ) );
      end
    end
  end
  
  function [vPath vLogLik] = viterbihelper(condLogLik,logStartP,logTransP)
    % Run the viterbi algorithm / max-product given probabilities
    % [vPath vLogLik] = viterbihelper(condLogLik,logTransP,logStartP)
    
    if exist('fastviterbi','file')==3 % Use the mex file (mex/dll==3)
      [vPath vLogLik] = fastviterbi(condLogLik,logStartP,logTransP);
    else
      [numHids numTimeSteps] = size(condLogLik);
      paths = zeros(numHids,numTimeSteps);
      
      paths(:,1) = 1:numHids;
      message = logStartP(:) + condLogLik(:,1);
      M = bsxfun(@plus,logTransP,message);
      [maxM paths(:,1)] = max(M,[],1); % Max over previous states
      for t=2:numTimeSteps
        message = condLogLik(:,t) + maxM(:);
        M = bsxfun(@plus,logTransP,message);
        [maxM paths(:,t)] = max(M,[],1); % Max over previous statess
      end
      [vLogLik maxPathInd] = max(message);
      vPath = paths(maxPathInd,:);
    end
  end
  function hidPath = samplepathhelper(numTimeSteps,sP,tP)
    % Sample a state sequence given start and transition probabilities
    % hidPath = HMM.samplepathhelper(numTimeSteps,sP,tP)
    
    if exist('fastsamplepath','file')==3 % Use the mex file (mex/dll==3)
      unifHelperSample = rand(1,numTimeSteps);
      hidPath = fastsamplepath(unifHelperSample,sP(:),tP');
    else
      assert(size(sP,1)==1 && ndims(sP)==2,'sP should be 1xnumHids');
      numH = length(sP);
      assert(ndims(tP)==2 && all(size(tP)==numH),...
        'tP should be numHidsxnumHids');
      hidPath = zeros(1,numTimeSteps);
      hidPath(:,1) = find(multisample(sP(:),1));
      for t=2:numTimeSteps
        p = tP(hidPath(:,t-1),:); % Row
        hidPath(:,t) = find(multisample(p(:),1));
      end
    end
  end
end
end
