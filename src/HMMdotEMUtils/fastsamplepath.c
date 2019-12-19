/* Sample a sequence of hidden states represented as a logical matrix
 * hidPath = fastsamplepath(unifHelperSample,sPT,tPT);
 * unifHelperSample is 1xnumTimeSteps uniform 0-1 numbers to avoid using stdlib
 * Expects transposes of transitionP for speed.
 */

#include "mex.h"
#include "matrix.h"

/* Generate a random index given a column probability distribution */
mwSize fsp_randind(double * p,mwSize numHids,double r)
{
  mwSize h,i=0;
  double sum=0;
  /* Previously uniform sampling was done here internally:
   * double r=rand()/(RAND_MAX+1.0) */
  
  for(h=0;h<numHids;h++){
    sum+=p[h];
    if(sum>=r){
      break;
    }else{
      i++;
    }
  }
/*
  mexPrintf("r=%f,i=%d\n",r,i);
*/
  return i;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *unifHelperSample,*sPT,*tPT; /* Expects transposes for speed */
  mwSize numTimeSteps,numHids;
  double *hidPath;
  mwSize t,h,i;
  
  /* Check inputs */
  if (nrhs!=3) mexErrMsgTxt("Expects exactly 3 inputs");
  if (nlhs>1) mexErrMsgTxt("Too many outputs");
  for(i=0;i<nrhs;i++){
    if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) )
      mexErrMsgTxt("All inputs must be full,real,double");
  }
  
  numTimeSteps = mxGetN(prhs[0]);
  if(mxGetM(prhs[0])!=1){
    mexErrMsgTxt("Helper uniform samples must be row vector");
  }
  unifHelperSample = mxGetPr(prhs[0]);
  if(numTimeSteps<0) mexErrMsgTxt("numTimeSteps must be positive");
  numHids = mxGetM(prhs[1]);
  if(mxGetN(prhs[1])!=1)
    mexErrMsgTxt("Expects startP and transitionP to be transposed!!!");
  if(mxGetM(prhs[2])!=numHids || mxGetN(prhs[2])!=numHids) 
    mexErrMsgTxt("Sizes of startP and transitionP not compatible");
  
  /* Create output */
  plhs[0] = mxCreateDoubleMatrix(1,numTimeSteps,mxREAL);
  if (numTimeSteps<=0) return;
  /* Get pointers */
  
  sPT = mxGetPr(prhs[1]);
  tPT = mxGetPr(prhs[2]);
  hidPath = mxGetPr(plhs[0]);
  
  /* Loop and sample */
  i = fsp_randind(sPT,numHids,unifHelperSample[0]);
  hidPath[0] = (double) i+1;
  for(t=1;t<numTimeSteps;t++){
    i = fsp_randind(tPT+i*numHids,numHids,unifHelperSample[t]);
    hidPath[t] = (double) i+1;
  }
  
}
