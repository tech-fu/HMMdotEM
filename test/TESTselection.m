function TESTselection(DATADIM,trueNumHids) % Needs to be a function!
% Test the model selection capabilities
if nargin<2
  trueNumHids = 2;
end
if nargin<1
  DATADIM=2;
end

trueGM = GaussMixHMM.randgaussmix(DATADIM,trueNumHids);
trueGM.assumeEquilibrium = true; % Make sure samples are from equilibrium

numDataSegs = 3;
maxData = 10^4; % Can jack up to 1 million and is still fast
minData = round(.9*maxData);
dataSegs = cell(numDataSegs,1);
truePaths = cell(numDataSegs,1);
randi = @(x) ceil(rand()*x);
for seg = 1:numDataSegs
   [dataSegs{seg} truePaths{seg}] =...
                            trueGM.sampledata(minData+randi(maxData-minData));
end

testData = trueGM.sampledata(maxData);

%% Try various models
minHids = max(1,trueNumHids - 1 - ceil(trueNumHids/2));
maxHids = trueNumHids+1+ceil(trueNumHids/2);
bic = zeros(maxHids,1);
icl = zeros(maxHids,1);
testLogLik = zeros(maxHids,1);
for h=minHids:maxHids
  GM = GaussMixHMM(dataSegs{1},h);
  GM = GM.constructparamgroup(dataSegs);
  dummyFM = FunMinimizer('fminsearch');
  EM = HMMEM(GM,dummyFM);
  
  EM.testall;
  EM.displayForwardBack = false;
  EM.displayHids = false;
  EM.displayData = true;
  EM.displayViterbi = false;
  EM.displayModel = false;
%   EM.displaySteps = false;
%   EM.displayAny = false; % Uncomment to see things
  EM.verbose = false; % Printing slows it down
  EM.maxSteps = 50;
  EM.logLTol = 10^-5;
  EM.run;
  
  if trueNumHids<=3 && maxData<=100 && h==trueNumHids 
    vP = [EM.Model(1).computeviterbi.viterbiPath];
    figure(12);clf;hold on;
    title('Path match? (Note: may be permuted)');
    plot(vP,'r-o');
    plot(truePaths{1},'k*');
    set(gca,'YLim',[.5 h+.5],'XLim',[-1 length(vP)+1]);
    legend('true','predicted');
    xlabel('time');
    ylabel('hidden state');
  end
  for seg=1:EM.numModels
    bic(h) = bic(h) + EM.Model(seg).bic;
    icl(h) = icl(h) + EM.Model(seg).icl;
  end
  testLogLik(h) = EM.Model(1).testloglik(testData);
  
  plotSelection(h);
  drawnow;
end

%% Plot selection criteria
function plotSelection(nH)
  figure(11);clf;hold on;
  plot(1:nH,bic(1:nH),'bo-');
  plot(1:nH,icl(1:nH),'r*-');
  plot(1:nH,testLogLik(1:nH),'k+-');
  title(['Model Selection for dataDim=' num2str(DATADIM)]);
  xlabel('number of hidden');
  ylabel('score');
  set(gca,'XTick',1:h);
  legend({'bic','icl','test data'},'Location','Best');
end
end
