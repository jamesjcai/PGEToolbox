#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"

#define FALSE  0
#define TRUE   1

/* compute the convolution of two sequences or filter one from the other */
void convolute(double *in1, int len1,double *in2,int len2, double *out)
{
    int i,j,k,m;
    double acc;
    
   m = 0;
   for (i=0; i<len1+len2-1; i++) {
        acc = 0;
        if (i > len1 -1) {
            j = len1 -1;
            m++;
        }
        else j = i;
        k = m;
        for (; j>=0 ; j--) {
            acc += in1[j] * in2[k];
            k++;
            if (k > len2 -1) break;
        }
        out[i] = acc;
    }        
}

void stirling(int n, double *out)
{
	int i,j, len1=2, len2=2;
	double in1[2]={1.0,-1.0};
	double in2[2];
	double *temp, *temp1;

	/* temp initialization */
	temp = malloc(len1*sizeof(double));
	for (j=0;j<len1;j++) temp[j] = in1[j];

	for (i=2;i<=n;i++) {
		temp1 = malloc((len1+1)*sizeof(double));
		in2[0] = 1.0;
		in2[1] = (float) -i*1.0;
		convolute(temp,len1,in2,len2,temp1);
		/* exchange temp to temp1 */
		len1 = len1 + 1;
		/* resize temp1 */
		free(temp);
		temp = malloc(len1*sizeof(double));
		for(j=0;j<len1;j++) temp[j] = temp1[j];
	}
	for(j=0;j<=n;j++) out[j]=temp1[j];
	free(temp1);
	free(temp);
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int n;
  double *x;
    
  if(nrhs!=1)
  mexErrMsgTxt("one input parameter is required.");
  
  if(nlhs!=1) 
  mexErrMsgTxt("One output required.");
    
  /* get the dimension of the data */
  n = (int) mxGetScalar(prhs[0]);
  
  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(1,n+1,mxREAL);
  x = mxGetPr(plhs[0]);
  
  /*  call the C subroutine */
  stirling(n,x);
}
