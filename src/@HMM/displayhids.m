function displayhids(obj,titlePrefix,useViterbi)
    % Plot the the most probable state (and sequence) using colors
    % obj.plotdata(titlePrefix,useViterbi)
    %
    % titlePrefix is a string that will be prepended to a title
    % useViterbi is true (Default) if you want to see most probable path
    %            the path will need to be computed
    assert(isscalar(obj),'cannot display hids for an array of HMMs');
    if nargin<2
      titlePrefix = '';
    end
    if nargin<3
      useViterbi = true;
    end
%     T = min(obj.numTimeSteps,obj.MAXDISPLAYTIME);
    T = obj.numTimeSteps;
    figtitle = @(s) title([titlePrefix s]);
    %% Prepare
    obj = obj.computeresp;
    if useViterbi
      zI = obj.computeviterbi.viterbiPath;
    end
    
    K = obj.numHids;
    allMarkers = allmarkers;
    numMarkers = size(allMarkers,1);
    needMarkers = K-numMarkers;
    if needMarkers>0 % If the markers are not enough pad with one value
      allMarkers(numMarkers:numMarkers+needMarkers,:) = 'k.';
    end
    markSize = obj.displaymarkersize;
    
    %% Plot resp
    subplot(1+double(useViterbi),1,1);hold on;
    for k=1:K
      imagesc(obj.computeresp.resp,[0 1]);colormap(gray);
    end
    xlabel('time'); ylabel('state');
    if obj.numHids~=1
      set(gca,'XLim',[1 T],'YLim',[1 obj.numHids]);
    end
    axis ij;
    figtitle('P(Z_t | data)');
    
    if useViterbi
      %% Show the viterbi state path
      subplot(2,1,2);hold on;
      plot(zI,'k.-')
      for k=1:K
        isK = zI==k;
        tK = find(isK);
        plot(tK,zI(isK),allMarkers(k,:),'MarkerSize',markSize);
      end
      xlabel('time'); ylabel('state');
  %     legend('Path');
      set(gca,'XLim',[1 T],'YLim',[.5 K+.5]);
      axis ij
      title(sprintf('Most probable path\nlogL = %.4e',obj.viterbiLogLik));
    end
  end % END OF method displaydata
