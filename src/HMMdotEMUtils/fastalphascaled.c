/* Compute the scaled forward quantities
 * This is intended for speedup.
 */

#include "mex.h"

void fas_todistribution(double * alphaPtr,double * scalingPtr,mwSize numHids){
  /* alphaPtr[0]-[numHids-1] will be normalized
   * scalingPtr[0] will be the normalization constant
   */
  mwSize h;
  scalingPtr[0] = 0;
  for(h=0;h<numHids;h++) scalingPtr[0] += alphaPtr[h];
  for(h=0;h<numHids;h++) alphaPtr[h] = alphaPtr[h]/scalingPtr[0];
}

void fas_computetrans(double * trans,double * alphaS,double * tP,mwSize numHids){
  /* tP'*alphaS */
  mwSize hr,hc,hInd;
  for(hc=0;hc<numHids;hc++){
    trans[hc] = 0;
    hInd = hc*numHids;
    for(hr=0;hr<numHids;hr++){
      trans[hc] += alphaS[hr]*tP[hr+hInd];
    }
/*
    mexPrintf("tr[%d]=%f, ",hc,trans[hc]);
*/
  }
}
        
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwSize T,numHids,t,h,tInd,linInd;
  double *condLik,*sP,*tP,*alphaS,*scaling;
  double sum,*trans;
  int i;
  
  /* Check inputs */
  if(nrhs!=3) mexErrMsgTxt("Expects exactly 3 inputs");
  if(nlhs>2) mexErrMsgTxt("Too many outputs");
  for(i=0;i<nrhs;i++){
    if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) )
      mexErrMsgTxt("All inputs must be full,real,double");
  }
  
  numHids = mxGetM(prhs[0]);
  if(mxGetM(prhs[1])!=1 || numHids!=mxGetN(prhs[1]) )
    mexErrMsgTxt("size(sP)~=[numHids 1]");
  if(numHids!=mxGetM(prhs[2]) || numHids!=mxGetN(prhs[2]))
    mexErrMsgTxt("size(tP)~=[numHids numHids]");
  
  T = mxGetN(prhs[0]);
  
  /* Create outputs */
  plhs[0] = mxCreateDoubleMatrix(numHids,T,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1,T,mxREAL);
  
  if (T==0 || numHids==0) return;
  
  /* Get pointers */
  condLik = mxGetPr(prhs[0]);
  sP = mxGetPr(prhs[1]);
  tP = mxGetPr(prhs[2]);
  alphaS = mxGetPr(plhs[0]);
  scaling = mxGetPr(plhs[1]);
  
  trans = mxMalloc(numHids*sizeof(double));
          
  /* Initial step */
  for(h=0;h<numHids;h++){
    alphaS[h] = condLik[h]*sP[h];
/*
    mexPrintf("%f, ",alphaS[h]);
*/
  }
  fas_todistribution(alphaS,scaling,numHids);
  fas_computetrans(trans,alphaS,tP,numHids);
/*
  mexPrintf("%f %f\n",alphaS[0],alphaS[1]);
*/
  tInd = 0;
  for(t=1;t<T;t++){
    tInd = t*numHids;
/*
    mexPrintf("tInd=%d, ",tInd);
*/
    
    for(h=0;h<numHids;h++){
      linInd = h+tInd;
      alphaS[linInd] = condLik[linInd]*trans[h];
    }
    fas_todistribution(alphaS+tInd,scaling+t,numHids);
    fas_computetrans(trans,alphaS+tInd,tP,numHids);
/*
    mexPrintf("\n");
*/
  }
}
