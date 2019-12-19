classdef HMMEM < handle 
% Handles EM computations for an HMM based on a HMM model.
% Stores likelihoods for each EM step
% It also has a property that says which function to use to optimize the M-Step
% and stores the gain from applying this minimizer to the data term
% (The transition and initial terms have closed-form solutions)
%
% The class inherits from handle, because this way you can kill a process
% and rerun the algorithm from where you left off.
% However being a handle makes this like a pointer, so 
% One = HMMEM;
% Two = One;
% One.Model = HMM(...); % Will also change the Model property of Two;
%
% This class can also take in HMM arrays
% they must be "param group" arrays
% that pass the HMM.isparamgroup test.
% Also see the HMM.constructparamsgroup static method

%% Constants
properties(Constant)
  % Temporary figures that try not to interfere with existing open ones.
  % See method genfignum
  FIGNUMS = 301:700;
  
  % Data will be plotted on figures > DATAFIG
  % For params groups all plots >DATAFIG might be used
  % For one Model only DATAFIG figure will be used.
  DATAFIG = 701;
  
  % Hidden state inference will be plotted on figures > HIDSFIG
  % Acts analogously to DATAFIG (see that)
  HIDSFIG = 801;
  
  MODELPFIG = 101; % Figure for plotting model parameters
  STEPSFIG = 102; % Figure for progress of optimization
  FORWARDBACKFIG = 110; % Figure for plotting forward-back quantities
  
end

%% Main
properties(SetAccess=public)
  % The HMM instance on which this will be run
  % Can also be an array of HMM's that pass the test
  % assert(Model.isparamgroup)
  % indicating that they have shared parameters (and possibly different data)
  Model
  
  % A function that solves the mStep optimization for the data term
  % Full freedom. Must be callable by 
  % feval(MMinimizer,qFun,condParams)
  % i.e. can optimize a function with a matrix argument
  % The variables qFun,condParams
  % are quantities computed internally
  %
  % Use FunMinizer class if you want to control optimization parameters
  % example:
  %   >> MMinimizer = FunMinimizer('fminsearch');
  %   >> MMinimizer.funParams.Display = 'off';
  MMinimizer;
end

%% Internal variables
properties(SetAccess=private)
  
  %% Keep track of steps
  
  % Number of EM steps that have been executed
  % Will increment when method "run" is called
  % You need to increment this in order to define a new em step
  % otherwise any variables that contain the progress from 
  % the previous step will be overwritten.
  numSteps = 0;
  
  % Log-L for each member of the param group at each EM step
  % paramGroupLogLik(emStep,d) is Model(d).hmmLogLik at emStep
  paramGroupLogLik 
  
  % Values of the function q (data term only) that is computed in the E Step
  % Two values per M-Step the first is the value before optimization
  % the second is after optimization
  % qDataTerm*[1;-1] will give the increment at each step resulting
  % from applying MMinimizer
  qDataTerm
  
end

properties(Dependent)
  % Log-Likelihoods, one for each EM-Step. On all the data
  % Computed as sum(obj.paramGroupLogLik,2)
  modelLogLik
  
  hasModelParamGroup % Whether Model property is a param group (must be true)
  numModels % number of Models (the number of HMM's in the Model property)
  
  isConverged % Whether the convergence conditions have been met
end

%% Parameters that can be altered directly
properties
  %% Flags
  
  verbose = true; % Print stuff?
  displayAny = true; % Display stuff? (e.g. Figures)
  
  % Whether to do checks to make sure algorithm is running properly
  % if true will run assertion statements that might require significant
  % computation, and will not be necessary if code is ok.
  testing = false;
  
  % Theses can be used to control plotting in more detail
  % Setting property display to false disables all plotting
  % but does not effect these values.
  
  displaySteps = true; % Display log-L and minimizer gain for all em steps?
  displayData = true; % Display data?
  displayHids = true; % Display hidden states?
  displayViterbi = true; % Display most probable path? Will decrease speed.
  displayModel = true; % Display the model parameters?
  
  displayForwardBack = false; % Display quantities in forward-back computation
  
  %% Figures
  
  % Helps to open new figures and keep track of open ones
  currentFigNum = HMMEM.FIGNUMS(1);
  
  %% Convergence
  
  logLTol = 10^-10; % Tolerance for convergence of EM.
  maxSteps = 100; % Maximum number of EM steps to run.
end

%% CONSTRUCTOR AND INITIALIZATION
methods
  function obj = HMMEM(Model,MMinimizer)
    % Construct an HMMEM based on a given model and MMinimizer
    % obj = HMMEM(Model,MMinimizer)
    
    if nargin>0
      obj.Model = Model;
      obj.MMinimizer = MMinimizer;
    end
  end
end

%% GET/SET
methods
  %% Set
  function set.currentFigNum(obj,figNum)
    % Sets the current figure number and makes sure it is one of FIGNUMS
    if ismember(figNum,HMMEM.FIGNUMS)
      obj.currentFigNum = figNum;
    end
  end
  function set.Model(obj,Model)
    assert(isa(Model,'HMM'),'Cannot set Model property to non-hmm class');
    if obj.testing %#ok<MCSUP> OK since "if []" will be ok
      assert(Model.isparamgroup,'Model is not a param group type HMM')
    end
    obj.Model = Model;
  end
  function set.MMinimizer(obj,MMinimizer)
    assert(canfeval(MMinimizer),...
        'MMinimizer must be passable to feval');
    obj.MMinimizer = MMinimizer;
  end
  function testall(obj)
    % Set obj.test = true and obj.Model.test = true
    obj.testing = true;
    [obj.Model.testing] = deal(true);
  end
  
  %% Get
  function bool = get.hasModelParamGroup(obj)
    bool = obj.Model.isparamgroup;
  end
  
  function numModels = get.numModels(obj)
    numModels = length(obj.Model);
  end
  
  function modelLogLik = get.modelLogLik(obj)
    modelLogLik = sum(obj.paramGroupLogLik,2);
  end
  
  function bool = belowtol(obj)
    % Check whether the last change in modelLogLik was below logLTol
    % bool = obj.belowtol()
    L = obj.modelLogLik;
    if ~isempty(L) && L(end)~=0
      bool = abs(L(end)-L(end-1))/abs(L(end-1)) < obj.logLTol;
    end
  end
  function bool = maxedsteps(obj)
    % Whether numSteps >= maxSteps
    % bool = obj.maxedsteps();
    bool = obj.numSteps>=obj.maxSteps;
  end
  function bool = get.isConverged(obj)
    bool = obj.numSteps>1 && (obj.maxedsteps || obj.belowtol);
  end
end

%% EM
methods
  function estep(obj)
    % Complete a full E-step. Updates "modelLogLik" and "paramGroupLogLik"
    % obj.estep()
    %
    % Updates Model "resp" and "transitionStats" needed in M-step
    
    numModels = obj.numModels;
    
    for i=1:numModels
      obj.Model(i) = obj.Model(i).forwardback;
      obj.Model(i) = obj.Model(i).computehmmloglik;
      obj.Model(i) = obj.Model(i).computeresp;
      obj.Model(i) = obj.Model(i).computetransitionstats;
      obj.paramGroupLogLik(obj.numSteps,i) = obj.Model(i).hmmLogLik;
      if obj.displayAny
        obj.displayforwardback;
      end
      % Free memory from varaibles unneeded in M step
      obj.Model(i) = obj.Model(i).emptyforwardback;
    end
    modelLogLik = obj.modelLogLik;
    
    if obj.numSteps>1 
      before = modelLogLik(obj.numSteps-1,1);
    else
      before = -inf;
    end
    current = modelLogLik(obj.numSteps,1);
    obj.fprintf('| logL = %+-.4e , dL = %+-.4e |',current,current-before);
    
  end
  
  function mstep(obj)
    % Runs an M-Step. Changes Model.condParams and updates property "qDataTerm"
    % obj.mstep()
    % Must have done E-Step first
        
    % Some simplfications of code
    numHids = obj.Model(1).numHids;
    numModels = obj.numModels;
    
    % Update the transition probs first to clear all the unneeded quantities
    % such as alphaS and betaS (modifying the model empties them)
    transitionStats = reshape([obj.Model.transitionStats],...
                                      numHids,numHids,numModels);
    transitionStats = sum(transitionStats,3);
    respOne = zeros(1,numHids);
    for i=1:numModels
      respOne = respOne + obj.Model(i).resp(:,1)';
    end
    sP = todistribution(respOne,'row');
    tP = todistribution(transitionStats,'row');
    % Need to modify only the first one, since Model.Params is a shared handle
    
    % Compute the data term update
    [newParams oldF newF] = obj.Model.fitcondparams(obj.MMinimizer);
    
    % Make the updates
    obj.Model = obj.Model.setstartp(sP);
    obj.Model = obj.Model.settransitionp(tP);
    if obj.testing
      tempL = sum(obj.Model.paramgrouploglik);
      tempDiff = tempL-obj.modelLogLik(end);
      obj.fprintf('| sP&tP %+-.4e |',tempDiff);
% Sometimes goes down by really tiny. assert(tempDiff>=0,'Log-L going down')
    end
    
    obj.Model(1) = obj.Model(1).setcondparams(newParams);
    
    % Progress
    obj.qDataTerm(obj.numSteps,:) = [oldF newF];
    
    % Verbose
    obj.fprintf('| MMinimizer = %+-.4e |',oldF-newF);
  end

  function run(obj)
  	% Run EM till convergence
  	% obj.run()
  	% Will increment the numSteps property.
    % See explanation about incrementing numSteps in the
    % numSteps property comments
    while ~obj.isConverged
      obj.numSteps = obj.numSteps+1;
      obj.fprintf('EM Step %d |',obj.numSteps);
      
      % E + logLik
      obj.estep;
      if obj.displayAny % This has to be before M-step
        obj.displaydata;
        obj.displayhids;
      end
      
      % M + functionValues
      obj.mstep;
      
      if obj.displayAny
        obj.displaysteps;
        obj.displaymodel;
      end
      
      obj.fprintf('\n'); % All statements inside should have been without \n
      if obj.displayAny
        drawnow;
      end
    end
  end
end

%% DISPLAY
methods
  function figNum = genfignum(obj)
    % Return the current figure number and increment it
    % so that the next call will get an unused one
    % obj.genfignum()
    figNum = obj.currentFigNum;
    FN = HMMEM.FIGNUMS; % Windows issue: HMMEM.FIGNUMS(end) doesn't work
    if obj.currentFigNum==FN(end)
      obj.currentFigNum = FN(1);
    else
      obj.currentFigNum = obj.currentFigNum+1;
    end
  end
  
  function s = figtitle(obj,titleStr)
    % Call title but prepend a string that tells about the current step
    % figtitle(obj,titleStr)
    % calls "title(s)", where s = [<EM info> titleStr]
    %
    % s = figtitle(obj,titleStr)
    % Doesn't call title, but returns the string so it can be altered
    
    s = sprintf('EM Step %d\n%s',obj.numSteps,titleStr);
    if nargout==0
      title(s);
    end
  end
  
  function displayforwardback(obj)
    % Display forward-back algorithm quantities
    % obj.displayforwardback()
    if obj.displayForwardBack && isscalar(obj.Model) % Cannot do it for non-scalar
      
      figure(obj.FORWARDBACKFIG);clf;
      subplot(2,2,1);
      imagesc(obj.Model.condLik,[0 1]);colormap(gray);
      xlabel('time');ylabel('state');
      obj.figtitle('condLik = P(data_t | Z_t)');
      
      subplot(2,2,2);
      imagesc(obj.Model.alphaS,[0 1]); colormap(gray);
      xlabel('time');ylabel('state');
      obj.figtitle('\alpha_{scaled} = P(Z_t | data_{1:t})');

      subplot(2,2,3);
      bar(obj.Model.scaling); colormap(gray);
      xlabel('time');ylabel('Prob');
      obj.figtitle('scaling = P(data_t | data_{1:t-1})');

      minBS = min(obj.Model.betaS(:));
      maxBS = max(obj.Model.betaS(:));
      subplot(2,2,4);
      imagesc(obj.Model.betaS); colormap(gray);
      xlabel('time');ylabel('state');
      NEWLINEASCII = 10;
      obj.figtitle(...
       ['\beta_{scaled} = P(data_{t+1:T} | Z_t)/P(data_{t+1:T} | data_{1:t})'...
        NEWLINEASCII 'min: ' num2str(minBS) ' max: ' num2str(maxBS)]);
    end
  end
  
  function displaysteps(obj)
    % Display progress of EM over the number of steps performed
    % obj.displaysteps()
    if obj.displaySteps
      figure(obj.STEPSFIG);clf;
      % Show the modelLogLik and qDataTerm values
      subplot(1,2,1);
      plot(obj.modelLogLik,'.-');
      xlabel('EM step');
      ylabel('log-likelihood');
      obj.figtitle('Log-likelihood of model');

      subplot(1,2,2);
      bar(obj.qDataTerm*[1;-1]);
      xlabel('EM step'); ylabel('Difference in data term');
      obj.figtitle('Gain from condParam update');
    end
  end
  function displaydata(obj)
    % Display the data and show the most probable state using colors
    % obj.displaydata()
    % opens an additional figure for each member of the param group
    if obj.displayData
      for i=1:obj.numModels
        figure(obj.DATAFIG+i-1);clf;
        obj.Model(i).displaydata( ...
          obj.figtitle(['Model ' num2str(i) 10]) , obj.displayViterbi);
      end
    end
  end
  function displayhids(obj)
    % Display the hidden states and show the most probable state using colors
    % obj.displayhids()
    % opens an additional figure for each member of the param group
    if obj.displayHids
      for i=1:obj.numModels
        figure(obj.HIDSFIG+i-1);clf;
        obj.Model(i).displayhids(...
          obj.figtitle(['Model ' num2str(i) 10]), obj.displayViterbi );
      end
    end
  end
  
  function displaymodel(obj)
    % Display the Model parameters
    % obj.displaymodel()
    if obj.displayModel
      figure(obj.MODELPFIG);clf;
      obj.Model(1).displayhmm;
      subplot(1,3,2);hold on;
      obj.figtitle('Transition Probs');
    end
  end
  
end

%% FLAGS
methods
  function fprintf(obj,varargin)
    % Call fprintf if obj.verbose
    % obj.fprintf(varargin)
    if obj.verbose
      fprintf(varargin{:});
    end
  end
end
end
