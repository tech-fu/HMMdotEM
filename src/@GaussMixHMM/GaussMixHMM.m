classdef GaussMixHMM < HMM
% Gaussian Mixture HMM, inherits from HMM
% Mainly used for testing HMM related function
% has simple methods

properties(Dependent)
  means
  covs
  
end

%% Constructor
methods
  function obj = GaussMixHMM(data,numHids,condArg)
    
    if nargin==0
      HMMArgs = {};
    else
      dataDim = size(data,1); 
      if nargin<3 || isempty(condArg)
        condArg = dataDim*(dataDim+3)/2;
      end
      HMMArgs = {data,numHids,condArg};
    end
    obj = obj@HMM(HMMArgs{:});
    
  end
  
end

%% Gaussian-specific
methods
  function logL = condloglikfun(obj,params,data)
    % Return the log-likelihood under a 2D Gaussian model
    % logL = gaussiancondloglik(params,data)
    % params(:,i) - [m_1 ; m_2 ; s_11 ; s_21 ; s_22] == [m(:);s([1 2 4])]; (mean;var)
    %          Note only one of the antidiagonal entries is given
    % data - 2 x numPoints
    %
    % logL - 
    % if params has many columns it will return a larger array with a row
    % of log-likelihoods for every column of params

    numClusters = size(params,2);
    [dataDim numPoints] = size(data);

    logL = zeros(numClusters,numPoints);
    meanInd = 1:dataDim;
    isUpper = triu(true(dataDim));
    isAbove = triu(true(dataDim),1);
    for i = 1:numClusters
      m = params(meanInd,i);
      V = ones(dataDim);
      V(isUpper) = params(meanInd(end)+1:end,i);
      VT = V';
      V(~isUpper) = VT(~isUpper);
      u = bsxfun(@minus,data,m);
      term = V\u;
      dV = det(V);
      if dV<=0 % Optimizers will not know about restrictions
        logL(i,:) = -inf;
      else
        logL(i,:) = -reallog(2*pi) - 0.5*(  reallog(dV) + sum(u.*term,1)  );
      end
    end

  end
  function obj = initcondparams(obj,numCondParams) %#ok<INUSD>
    % Estimate/initialize conditional parameters
    % uses kmeans
    % obj = obj.initcondparams()
    
    assert(~isempty(obj.data),'cannot initialize with empty data');
    
    MAXT = 10^6;
    X = obj.data(:,1:min(obj.numTimeSteps,MAXT))';
    [ids m] = kmeans(X,obj.numHids);
    
    V = cov(obj.data',1);
    isUpper = triu(true(obj.dataDim));
    VP = V(isUpper)*ones(1,obj.numHids);
    cP = [m' ; VP];
    obj = obj.setcondparams(cP);
  end
  
  function [newParams oldVal newVal] = fitcondparams(obj,varargin)
    % Compute an analytic parameter update
    % [newParams oldVal newVal] = obj.fitcondparams(varargin)
    % 
    % Works with no arguments
    % the "varargin" is only to satisfy the superclass prototype
    
    assert(obj.isparamgroup,'This HMM is not a paramgroup');
    weights = [obj.resp]; % Works only for row-time data
    weights = bsxfun(@rdivide,weights,sum(weights,2)); % Need weights
    meansInd = 1:obj(1).dataDim;
    
    dataAll = [obj.data];
    isUpper = triu(true(obj(1).dataDim));
    newParams = obj(1).condParams;
    newParams(meansInd,:) = dataAll*weights';
    
    for h=1:obj(1).numHids
      z = bsxfun(@minus,dataAll,newParams(meansInd,h));
      C = bsxfun(@times,z,weights(h,:))*z';
      newParams(meansInd(end)+1:end,h) = C(isUpper);
    end
    oldVal = 0;
    newVal = 0;
  end
end

%% GET/SET
methods
  function m = get.means(obj)
    m = obj.condParams(1:obj.dataDim,:);
  end
  function C = get.covs(obj)
    dDim = obj.dataDim;
    indMat = ones(dDim);
    TMAT = true(dDim);
    isUpper = triu(TMAT);
    isLower = isUpper';
    indMat(isUpper) = 1:(dDim*(dDim+1)/2);
    indMatT = indMat';
    indMat(isLower) = indMatT(isLower);
    C = reshape(obj.condParams(dDim+indMat,:),dDim,dDim,[]);
  end
  
  function obj = namecondparams(obj)
    totalParams = obj.numCondParams;
    condParamNames = cell(totalParams,1);
    dDim = obj.dataDim;
    meanInd = 1:dDim;
    for i=meanInd
      condParamNames{i} = ['m' num2str(i)];
    end
    [v1 v2] = find(triu(true(dDim)));
    for i=dDim+1:totalParams
      condParamNames{i} =...
        ['v' num2str(v1(i-dDim)) num2str(v2(i-dDim))];
    end
    obj = namecondparams@HMM(obj,condParamNames);
  end
end

%% Sampling
methods
  function [data hidPath] = sampledata(obj,numData)
    % Sample some data from the model
    % [data hidPath] = obj.sampledata(numData)
    %
    % numData - number of data points
    % data - the resulting data
    % binHidMat - binary hidden path numHids x numData indicating the states
    
    dDim = obj.dataDim;
    numH = obj.numHids;
    hidPath = obj.samplepath(numData);
    data = zeros(dDim,numData);
    isHid = accumarray(hidPath(:),1:numData,[numH 1],@(x) {x},{[]});
    numOccur = cellfun(@length,isHid);
    for k=1:obj.numHids
      m = obj.means(:,k);
      V = obj.covs(:,:,k);
      S = sqrtm(V);
      data(:,isHid{k}) = randgauss(m,S,numOccur(k));
    end
  end
  
end

%% Display
methods
  function ellipse = plotclusters(obj,k)
    % Get the k'th or all ellipses based on each cluster
    % ellipse = plotclusters(obj,k) - for the k'th
    % ellipse = plotclusters(obj) - for all
    %
    % ellipse is a (1 or numHids)x(Enough points) matrix
    % ellipse is based on unit circle multiplied by standard deviation matrix
    if nargin==1
      k = 1:obj.numHids;
    end
    if nargout>0
      assert(isscalar(k),...
        'if you ask the function to return an ellipse k must be a scalar')
    end
    theta = 0:.02:2*pi;
    circle = [cos(theta);sin(theta)];
    for i=k
      m = obj.means(:,i);
      V = obj.covs(:,:,i);
      S = sqrtm(V);
      ellipse = bsxfun(@plus,S*circle, m);
      if nargout==0
        plot(ellipse(1,:),ellipse(2,:),'g-','LineWidth',2);
      end
    end
  end
  
  function displaydata(obj,titlePrefix,useViterbi)
    % (Redefined) Display data as clusters
    % obj.displaydata(titlePrefix,useViterbi)
    assert(isscalar(obj),'cannot display data for an array of HMMs');
    if nargin<2
      titlePrefix = '';
    end
    if nargin<3
      useViterbi = true;
    end
    
    dDim = obj.dataDim;
    if dDim==2
      zI = cell(2,1);
      figtitle = @(s) title([titlePrefix s]);
      %% Prepare
      allMarkers = allmarkers;
      [unused zI{1}] = max(obj.computeresp.resp,[],1);
      clear unused
      if useViterbi
        zI{2} = obj.computeviterbi.viterbiPath;
      end

      x = obj.data(1,:);
      y = obj.data(2,:);
      K = obj.numHids;

      numMarkers = size(allMarkers,1);
      needMarkers = K-numMarkers;
      if needMarkers>0 % If the markers are not enough pad with one value
        allMarkers(numMarkers:numMarkers+needMarkers,:) = 'k.';
      end

      for i=1:(useViterbi+1)
        %% Plot x vs y
        subplot(useViterbi+1,1,i);hold on;
        if i==2
          plot(x,y,'k.-')
        end
        for k=1:K
          isK = zI{i}==k;
          plot(x(isK),y(isK),allMarkers(k,:));
        end
        xlabel('x'); ylabel('y');
        if i==2
          figtitle('Data - Viterbi');
        else
          obj.plotclusters;
          title('Data - Clusters');
        end
      end
    else % dDim~=2
      displaydata@HMM(obj,titlePrefix,useViterbi);
    end
  end
end

%% STATIC
methods(Static)
  function GM = randgaussmix(dataDim,numHids)
    % Construct a random GaussMixHMM with empty data
    % obj = GaussMixHMM.randgaussmix(dataDim,numHids)
    isUpper = triu(true(dataDim));

    condParams = zeros(dataDim+nnz(isUpper),numHids);
    meanStd = 10;
    for h=1:numHids
      condParams(1:dataDim,h) = meanStd*randn(dataDim,1);
      V = wishrnd(eye(dataDim),dataDim);
      condParams(dataDim+1:end,h) = V(isUpper);
    end

    rnddist = @(nR,nC) todistribution(gamrnd(10,.1,nR,nC),'row');
    sP = rnddist(1,numHids);
    tP = rnddist(numHids,numHids);

    GM = GaussMixHMM(zeros(dataDim,0),numHids,condParams);
    GM = GM.setstartp(sP);
    GM = GM.settransitionp(tP);
  end
end
end
