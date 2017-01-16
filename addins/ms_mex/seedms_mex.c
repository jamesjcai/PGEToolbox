#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mex.h" 

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {
  unsigned int iseed;
  if(nrhs>1|nlhs>0) {
	  mexErrMsgTxt("Usage: SEEDMS(12345)\n");
  }  
  
  if(nrhs==1) {
        if (mxIsEmpty(prhs[0])|| (!mxIsNumeric(prhs[0]))){
             mexErrMsgTxt("NSAM must be an integer.");
        }
     iseed = (int)mxGetScalar(prhs[0]);
  }  else {      
     /* Simple "srand()" seed: just use "time()" */
     iseed = (unsigned int)time(NULL);
  }
  srand (iseed);
  mexPrintf ("%d\n", iseed);
  return;
}
