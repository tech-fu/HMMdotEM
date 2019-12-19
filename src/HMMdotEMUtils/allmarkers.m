function specs = allmarkers()
% Return all combinations of standard markers
% in non-uniform way (neighboring markers will not be similar at all)
% 
% specs = allmarkers()
% Can plot (x,y) using i'th spec:
%   >>plot(x,y,specs(i,:))

LINESPECCOLORS = ['r' 'b' 'c' 'm' 'g' 'y']';
LINESPECMARKERS = ['+' 'o' '*' '.' 'x' 's' 'd' '^' 'v' '<' '>' 'p' 'h']';

numColors = length(LINESPECCOLORS);
numMarkers = length(LINESPECMARKERS);
numAll = numColors*numMarkers;
% specsInd = zeros(numAll,2);
specs = char(zeros(numAll,2));

%% Go down each diagonal
currentInd = 0;
for i=1:numMarkers
  for j=1:numColors
    currentInd = currentInd+1;
    mInd = i+j-1 - numMarkers*(i+j-1>numMarkers);
    cInd = j;
    specs(currentInd,:) = [LINESPECMARKERS(mInd) LINESPECCOLORS(cInd)];
  end
end

