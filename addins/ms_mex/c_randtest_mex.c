/*
 * rand: Generates 5 numbers using standard "srand()/rand()" function
 *
 * SAMPLE OUTPUT:
 *   rand[0]= 824522256
 *   rand[1]= 1360907941
 *   rand[2]= 1513675795
 *   rand[3]= 1046462087
 *   rand[4]= 253823980
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mex.h" 

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {
  
   /* Simple "srand()" seed: just use "time()" */
   //unsigned int iseed = (unsigned int)time(NULL);
   //srand (iseed);

  /* Now generate 5 pseudo-random numbers */
  int i;
  for (i=0; i<5; i++)
  {
    mexPrintf ("rand[%d]= %u\n",
      i, rand ());
  }
  return;
} 

