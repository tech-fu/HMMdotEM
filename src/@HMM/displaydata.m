function displaydata(obj,titlePrefix,useViterbi)
    % Plot the data and show the most probable state sequence using colors
    % obj.plotdata(titlePrefix)
    %
    % titlePrefix is a string that will be prepended to a title
    % useViterbi is true (Default) if you want to see most probable path
    %            the path will need to be computed
    
    if nargin<2
      titlePrefix = '';
    end
    if nargin<3
      useViterbi = true;
    end
    assert(isscalar(obj),'cannot display data for an array of HMMs');
    
    figtitle = @(s) title([titlePrefix s]);
    %% Prepare
    
    if useViterbi
      zI = obj.computeviterbi.viterbiPath;
    else
      [unused zI] = max(obj.computeresp.resp,[],1);
      clear unused;
    end
    
    K = obj.numHids;
    allMarkers = allmarkers;
    numMarkers = size(allMarkers,1);
    needMarkers = K-numMarkers;
    if needMarkers>0 % If the markers are not enough pad with one value
      allMarkers(numMarkers:numMarkers+needMarkers,:) = 'k.';
    end
    markSize = obj.displaymarkersize;
    
    %% Plot x vs time
    dDim = obj.dataDim;
    for dimInd = 1:dDim
      x = obj.data(dimInd,:);
      xStr = ['x_' num2str(dimInd)];
      subplot(dDim,1,dimInd);hold on;
      if useViterbi
        plot(x,'k.-')
      end
      for k=1:K
        isK = zI==k;
        tK = find(isK);
        plot(tK,x(isK),allMarkers(k,:),'MarkerSize',markSize);
      end
      
      xlabel('time'); ylabel(xStr);
  
      if dimInd==1
        if useViterbi
          figtitle('Viterbi');
        else
          figtitle('Max resp');
        end
      end
    end
    
  end % END OF method displaydata
