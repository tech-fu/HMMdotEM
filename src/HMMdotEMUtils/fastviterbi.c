/* Compute the most probable path and its likelihood
 * [vPath vLogLik] = fastviterbi(condLogLik,logTransP,logStartP)
 */

#include "mex.h"
#include "matrix.h"

mwSize fv_max(double *x,mwSize numX){
  /* Loop through and find the index of the maximum value */
  mwSize i,maxInd;
  double maxVal;
  maxInd=0;
  maxVal = x[0];
  for(i=0;i<numX;i++){
    if(x[i]>maxVal) {
      maxInd=i;
      maxVal = x[maxInd];
    }
  }
  return maxInd;
}

void fv_margmat(double *M,double *logTransP,double *messages,mwSize numHids){
  /* Compute the matrix of messages to the next time step */
  mwSize hr,hc,hInd;
  for(hc=0;hc<numHids;hc++){
    hInd = hc*numHids;
    for(hr=0;hr<numHids;hr++,hInd++){
      M[hInd] = logTransP[hInd] + messages[hr];
    }
  }
}

void fv_maxmarg(double *maxM,mwSize *maxInd,double *M,mwSize numHids){
  /* Max-marginals */
  mwSize hr,hc,hInd;
  for(hc=0;hc<numHids;hc++){
    maxM[hc] = 0;
    hInd = hc*numHids;
    for(hr=0;hr<numHids;hr++){
      maxInd[hc] = fv_max(M+hInd,numHids);
      maxM[hc] = M[maxInd[hc]+hInd];
    }
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *condLogLik,*logTransP,*logStartP,*vPath,*vLogLik;
  mwSize numHids,numTimeSteps;
  double *paths,*messages,*maxM,*M;
  mwSize t,tInd,h,*maxInd,maxPathInd;
  int i;
  
  /* Check args */
  if (nrhs!=3) mexErrMsgTxt("Expects exactly 3 inputs");
  if (nlhs>2) mexErrMsgTxt("Too many outputs");
  for(i=0;i<nrhs;i++){
    if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) )
      mexErrMsgTxt("All inputs must be full,real,double");
  }
  
  numHids = mxGetM(prhs[0]);
  numTimeSteps = mxGetN(prhs[0]);
  if(mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=numHids)
      mexErrMsgTxt("logStartP not numHidsxnumHids");
  if(mxGetM(prhs[2])!=numHids || mxGetN(prhs[2])!=numHids)
      mexErrMsgTxt("logTransP not 1xnumHids");
  
  /* Create outputs */
  plhs[0] = mxCreateDoubleMatrix(1,numTimeSteps,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
  
  if (numTimeSteps==0 || numHids==0) return;
  
  /* Get pointers */
  condLogLik = mxGetPr(prhs[0]);
  logTransP = mxGetPr(prhs[2]);
  logStartP = mxGetPr(prhs[1]);
  
  vPath = mxGetPr(plhs[0]);
  vLogLik = mxGetPr(plhs[1]);
  
  /* Pass messages */
  paths = mxMalloc(numHids*numTimeSteps*sizeof(double));
  messages = mxMalloc(numHids*sizeof(double));
  M = mxMalloc(numHids*numHids*sizeof(double));
  maxM = mxMalloc(numHids*sizeof(double));
  maxInd = mxMalloc(numHids*sizeof(mwSize));
  
  for(h=0;h<numHids;h++){
    maxM[h] = logStartP[h];
  }
  
  for(t=0;t<numTimeSteps;t++){
    
    tInd = t*numHids;
    for(h=0;h<numHids;h++){
      messages[h] = maxM[h] + condLogLik[h+tInd];
    }
    fv_margmat(M,logTransP,messages,numHids);
    fv_maxmarg(maxM,maxInd,M,numHids);
    for(h=0;h<numHids;h++){
      paths[h+tInd] = 1.0+(double) maxInd[h];
    }
  }
  
  maxPathInd = fv_max(messages,numHids);
  vLogLik[0] = messages[maxPathInd];
  for(t=0;t<numTimeSteps;t++){
    vPath[t] = paths[maxPathInd + t*numHids];
  }
  
}
