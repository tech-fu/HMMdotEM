/* Compute the scaled forward quantities
 * This is intended for speedup.
 */

#include "mex.h"

void fbs_stepback(double * betaS,double * trans, double * tP,mwSize numHids){
  /* tP'*alphaS */
  mwSize hr,hc;
  for(hr=0;hr<numHids;hr++){
    betaS[hr] = 0;
    for(hc=0;hc<numHids;hc++){
      betaS[hr] += tP[hr+hc*numHids]*trans[hc];
/*
      mexPrintf("b[%d]=%f=%f*%f, \n",hr,betaS[hr],tP[hr+hc*numHids],trans[hc]);
*/
    }
  }
}
        
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwSize T,numHids,t,h,tInd,linInd;
  double *betaS,*condLik,*alphaS,*scaling,*tP;
  double *trans;
  int i;
  
  /* Check inputs */
  if(nrhs!=4) mexErrMsgTxt("Expects exactly 3 inputs");
  if(nlhs>1) mexErrMsgTxt("Too many outputs");
  for(i=0;i<nrhs;i++){
    if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) )
      mexErrMsgTxt("All inputs must be full,real,double");
  }
  
  numHids = mxGetM(prhs[0]);
  if(numHids!=mxGetM(prhs[3]) || numHids!=mxGetN(prhs[3]))
    mexErrMsgTxt("size(tP)~=[numHids numHids]");
  
  T = mxGetN(prhs[0]);
  if(mxGetM(prhs[1])!=numHids || mxGetN(prhs[1])!=T)
    mexErrMsgTxt("size(alphaS)~=[numHids T]");
  if(mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=T)
    mexErrMsgTxt("size(scaling)~=[1 T]");
  
  /* Create outputs */
  plhs[0] = mxCreateDoubleMatrix(numHids,T,mxREAL);
  
  if (T==0 || numHids==0) return;
  
  /* Get pointers */
  condLik = mxGetPr(prhs[0]);
  alphaS = mxGetPr(prhs[1]);
  scaling = mxGetPr(prhs[2]);
  tP = mxGetPr(prhs[3]);
  betaS = mxGetPr(plhs[0]);
  
  trans = mxMalloc(numHids*sizeof(double));
  
  /* Initial step */
  tInd = (T-1)*numHids;
  for(h=tInd;h<T*numHids;h++){
    betaS[h] = 1;
  }
  
  for(t=T-1;t>0;t--){
    for(h=0;h<numHids;h++){
      linInd = h+tInd;
      trans[h] = betaS[linInd]*(condLik[linInd]/scaling[t]);
/*
      mexPrintf("lI=%d trans[%d]=%f=%f*%f/%f, ",linInd,h,trans[h],
              betaS[linInd],condLik[linInd],scaling[t]);
*/
    }
    tInd = tInd-numHids;
    fbs_stepback(betaS+tInd,trans,tP,numHids);
/*
    mexPrintf("\n");
*/
  }
}
