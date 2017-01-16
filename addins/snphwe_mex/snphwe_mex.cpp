#include <stdio.h>
#include <algorithm>
#include <bitset>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <cassert>
#include "mex.h"
static double SNPHWEP(int obs_hets, int obs_hom1, int obs_hom2);



void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double p;
  int nhom1,nhet,nhom2;

  if(nrhs!=3) {
   mexPrintf("This function implements an exact SNP test of Hardy-Weinberg\n");
   mexPrintf("Equilibrium as described in Wigginton, JE, Cutler, DJ, and\n");
   mexPrintf("Abecasis, GR (2005) A Note on Exact Tests of Hardy-Weinberg\n");
   mexPrintf("Equilibrium. American Journal of Human Genetics. 76(5): 887–893\n\n");
   mexErrMsgTxt("Three inputs required, e.g.,SNPHWE_MEX(381,1467,23)\n[p_value]=SNPHWE_MEX(nhom1,nhet,nhom2);");
  } else if(nlhs>1) {
    mexErrMsgTxt("Usage: [p_value]=SNPHWE_MEX(nhom1,nhet,nhom2);\nToo many output arguments");
  }
  /* NSAM */
  if (mxIsEmpty(prhs[0])|| (!mxIsNumeric(prhs[0]))){
	  mexErrMsgTxt("Usage: [p_value]=SNPHWE_MEX(nhom1,nhet,nhom2);\nNHOME1 must be an integer.");
  }
  nhom1 = (int)mxGetScalar(prhs[0]);
   /* NSAM */
  if (mxIsEmpty(prhs[1])|| (!mxIsNumeric(prhs[1]))){
    mexErrMsgTxt("Usage: [p_value]=SNPHWE_MEX(nhom1,nhet,nhom2);\nNHET must be an integer.");
  }
  nhet = (int)mxGetScalar(prhs[1]);
   /* NSAM */
  if (mxIsEmpty(prhs[2])|| (!mxIsNumeric(prhs[2]))){
    mexErrMsgTxt("Usage: [p_value]=SNPHWE_MEX(nhom1,nhet,nhom2);\nNHOME2 must be an integer.");
  }
  nhom2 = (int)mxGetScalar(prhs[2]);

  p=SNPHWEP(nhet, nhom1, nhom2);
  plhs[0] = mxCreateDoubleScalar(p);
  return;
}


/*
// This function implements an exact SNP test of Hardy-Weinberg
// Equilibrium as described in Wigginton, JE, Cutler, DJ, and
// Abecasis, GR (2005) A Note on Exact Tests of Hardy-Weinberg
// Equilibrium. American Journal of Human Genetics. 76(5): 887–893
//
// Written by Jan Wigginton
*/

double SNPHWEP(int obs_hets, int obs_hom1, int obs_hom2)
{
	if (obs_hom1 + obs_hom2 + obs_hets == 0 ) return 1;

	//if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0)
	//	error("Internal error: negative count in HWE test", 91);

	int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
	int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

	int rare_copies = 2 * obs_homr + obs_hets;
	int genotypes   = obs_hets + obs_homc + obs_homr;

	double * het_probs = (double *) malloc((size_t) (rare_copies + 1) * sizeof(double));
	//if (het_probs == NULL)
	//	error("Internal error: SNP-HWE: Unable to allocate array", 90);

	int i;
	for (i = 0; i <= rare_copies; i++)
		het_probs[i] = 0.0;

	/* start at midpoint */
	int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);

	/* check to ensure that midpoint and rare alleles have same parity */
	if ((rare_copies & 1) ^ (mid & 1))
	mid++;

	int curr_hets = mid;
	int curr_homr = (rare_copies - mid) / 2;
	int curr_homc = genotypes - curr_hets - curr_homr;

	het_probs[mid] = 1.0;
	double sum = het_probs[mid];
	for (curr_hets = mid; curr_hets > 1; curr_hets -= 2)
	{
	  het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
	/ (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
	  sum += het_probs[curr_hets - 2];

	  /* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
	  curr_homr++;
	  curr_homc++;
	}

	curr_hets = mid;
	curr_homr = (rare_copies - mid) / 2;
	curr_homc = genotypes - curr_hets - curr_homr;
	for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2)
	{
	  het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc	/((curr_hets + 2.0) * (curr_hets + 1.0));
	  sum += het_probs[curr_hets + 2];

	  /* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
	  curr_homr--;
	  curr_homc--;
	}

	for (i = 0; i <= rare_copies; i++)
	het_probs[i] /= sum;

	/* alternate p-value calculation for p_hi/p_lo
	double p_hi = het_probs[obs_hets];
	for (i = obs_hets + 1; i <= rare_copies; i++)
	 p_hi += het_probs[i];

	double p_lo = het_probs[obs_hets];
	for (i = obs_hets - 1; i >= 0; i--)
	  p_lo += het_probs[i];

	double p_hi_lo = p_hi < p_lo ? 2.0 * p_hi : 2.0 * p_lo;
	*/

	double p_hwe = 0.0;
	/*  p-value calculation for p_hwe  */
	for (i = 0; i <= rare_copies; i++)
	{
	  if (het_probs[i] > het_probs[obs_hets])
	continue;
	  p_hwe += het_probs[i];
	}

	p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

	free(het_probs);

	return p_hwe;
}



