function sample = multisample(probs,N)

% Sample from multiple multinomial distributions
% sample = multisample(probs,N(optional))
% Samples randomly from probs matrix
% returns sample - sparse
% Assumes columns are distributions
% N is the number of samples in a multinomial (default 1)
%if N>1 the function uses a for loop
%
% probs should be a proper distribution, but it is up to the user to check

[r c] = size(probs); 
probs = cumsum(probs,1);
probs(end,:) = 1; % Safety precaution
if nargin==1
  N=1;
end

if nargin==1 || N==1
  R =   rand(1,c);
  I = sum(bsxfun(@le,probs,R),1) + 1;
%     I = sum(probs<= ones(r,1)*rand(1,c),1)+1;
    sample = sparse(I,1:c,1,r,c);
else
  %Should call c function
    sample = sparse([],[],[],r,c,N*c);
    for i=1:N
      I = sum(bsxfun(@le,probs,rand(1,c)),1) + 1;
      sample = sample + sparse(I,1:c,1,r,c);
    end
end

