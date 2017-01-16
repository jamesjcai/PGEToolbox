#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "get_optval.h"
#include "hweStrata.h"


int main( int argc, char *argv[]) {


  /****Default settings for input options*****/
  int help=0;
  int ncolGeno;
  int maxstrata;
  int *genoRow;
  int maxcat;

  double pStop = 1.1;
  char *genoFileName = "geno.dat";
  char line[MAXLINE];

  FILE *genoFile;
  GENODAT genoDat;

  char * cmd = cmdString(argc, argv);

  if(get_optval(cmd,"h",OV_FLAG,&help) == 1)
	{
	  printHelp(argv[0], genoFileName, pStop);
	  exit(1);
	}


  if( get_optval(cmd, "geno",OV_STRING, &genoFileName) == -1)
	{
	  printf("Error: no value for option geno\n");
	  exit(1);
	}

  if( get_optval(cmd,"pstop",OV_DOUBLE, &pStop) == -1)
	{
	  printf("Error: no value for option pstop\n");
	  exit(1);
	}


  genoFile=fopen(genoFileName,"r");
  if (genoFile == NULL) {
	printf("Error: Geno file not found\n");
	exit(1);
  }


  ncolGeno = countColsFile(genoFile);

  /* subtract 1 for label in 1st col */

  ncolGeno = ncolGeno - 1;

  if( (ncolGeno % 3) != 0){
	printf("Error: Number genotype counts of geno data not divisible by 3\n");
	exit(1);
  }


  /* Assume that all rows have the same number of columns, so that
	 max strata is ncolGeno/3, for 3 SNP genotypes per stratum
  */


  maxstrata = ncolGeno / 3;
  genoDat.maxstrata = maxstrata;

  /* temp storage of genotypes */
  genoRow = (int *) malloc((maxstrata * 3) * sizeof(int));

  if (!genoRow){
    printf("Error: Allocation failure for genoRow\n");
    exit(1);
  }


  maxcat = 3*maxstrata;

  initGenoDat(&genoDat, pStop, maxstrata);

  printHeader(maxstrata);

  while(fgets(line, MAXLINE, genoFile)){


	getGenoDat(line, genoRow, maxcat, maxstrata, &genoDat);

	/* Test HWD Homogeneity Across Strata */

	hwdHomogExact(&genoDat);


	/* Exact Test HWE For All Strata */

	hweStratExact(&genoDat);

	
	/* Haldane's stratified test */

	/* hweStratHaldane(&genoDat); */


	/* Olson's  stratified test */

	/* hweStratOlson(&genoDat); */


	printHweStrat(&genoDat);
	fflush(stdout);

  }
  fclose(genoFile);
  return 1;
}

/***********************************************************************************/

void initGenoDat(GENODAT *genoDat, double pStop, int maxstrata){
  int j;

  genoDat-> pStop = pStop;


  genoDat->geno = imatrix(maxstrata, 3);

  genoDat->genoLabel = (char *) malloc(MAXLABELSIZE * sizeof(char));
  if(!genoDat->genoLabel){
    printf("Error: memory alloc failure for genoLabel");
    exit(1);
  }

  /* set values to defaults */

  genoDat-> pvalHomog = 1.0;
  genoDat-> phweCombined = 1.0;;
  genoDat-> stopPval = 0;
  genoDat-> pvalHaldane  = 1.0;
  genoDat-> zstatHaldane = 0.0;
  genoDat-> pvalOlson   = 1.0;
  genoDat-> zstatOlson  = 0.0;


  genoDat->phweStratum = (double *) malloc(maxstrata * sizeof(double));
  if(!genoDat->phweStratum){
    printf("Error: Allocation failure  for phweStratum\n");
    exit(1);
  }

  for(j=0; j<maxstrata; j++){
    genoDat->phweStratum[j] = 1.0;
  }

  genoDat->countA  =  (int *) malloc(maxstrata * sizeof(int));
  genoDat->countB  =  (int *) malloc(maxstrata * sizeof(int));
  genoDat->nSubj   =  (int *) malloc(maxstrata * sizeof(int));
  genoDat->hetLow  =  (int *) malloc(maxstrata * sizeof(int));
  genoDat->hetHi   =  (int *) malloc(maxstrata * sizeof(int));
  genoDat->h       =  (int *) malloc(maxstrata * sizeof(int));
  genoDat->hLow    =  (int *) malloc(maxstrata * sizeof(int));
  genoDat->hHi     =  (int *) malloc(maxstrata * sizeof(int));
  genoDat->csumh       =  (int *) malloc(maxstrata * sizeof(int));
  genoDat->rcsumHetHi  =  (int *) malloc(maxstrata * sizeof(int));
  genoDat->index       =  (int *) malloc(maxstrata * sizeof(int));
  genoDat->indexUpper  =  (int *) malloc(maxstrata * sizeof(int));

 
  genoDat->lnpObsVec   = (double *) malloc(maxstrata * sizeof(double));

 /* allocate pointers to rows */
  genoDat->lnprob = (double **) malloc(maxstrata * sizeof(double));       
  if (!genoDat->lnprob){
    printf("Error: Allocation failure 1 for lnprob\n");
    exit(1);
  }


  /* allocate rows and set pointers to them */
  for(j=0; j <maxstrata; j++){
    genoDat->lnprob[j] = (double *) malloc((MAXHETEROZYGOTES/2 + 1) * sizeof(double));
    if(!genoDat->lnprob[j]){
      printf("Error: Allocation failure 2 for lnprob\n");
      exit(1);
    }
  }

}


/***********************************************************************************/

void printHelp(char* argvFirst, char*genoFileName, double pStop){ 
    
  printf("\nCommand line options for %s, version %2.1f \n\n",argvFirst, VERSION);
  printf("  -h          : print this help menu\n\n");
  printf("  -geno <file>: input genotype file\n");
  printf("               (e.g., -geno geno.dat) default = %s\n\n", genoFileName);
  printf("  -pstop      : stop computing pval for HWE if p > pstop\n");
  printf("               (e.g., -pstop 0.5) default = %f\n", pStop);

}

/**************************************************************************/

int countColsFile(FILE *file){

  char *token;
  int ncol = 0;
  char line[MAXLINE];

  fgets(line,MAXLINE,file);

  token=strtok(line,SEPCHAR);
  while(token != NULL){
    ncol ++;
    token=strtok(NULL,SEPCHAR);
  }
  rewind(file);
  return(ncol);
}

/***********************************************************************************/


void getGenoDat(char *line, int *genoRow, int maxcat, int maxstrata, GENODAT *genoDat){

  char *token;
  int col;
  int indxStrata;
  int i,j,k,n;

  /* get geno label from 1st col */

  token = strtok(line,SEPCHAR);
    
  if(token !=NULL){
    strcpy(genoDat->genoLabel,  token);
  }
  
  col = 0;

  while(token != NULL){
      
    token=strtok(NULL,SEPCHAR);
      
    if(token != NULL){

      if(col >= maxcat ){
		printf("Error: Geno file has a row with > %d genotypes\n", maxcat);
		exit(1);
      }

      genoRow[col] = atoi(token);
      col ++;
    }
  }

  /* check genotype counts per stratum, to be sure total N > 0. If not
     then do not store genotype counts for a given stratum */
 
  indxStrata = 0;

  genoDat-> nstrata = 0;
 
  for(i=0; i < maxstrata; i++){
    n = 0;
    for(j=0; j < 3; j++){
      k = i*3 + j;
      n += genoRow[k];
    }
    
    /* temp check below to force 0 counts in geno */
    if(n >= 0){ 
      genoDat-> nstrata ++;
      for(j=0; j< 3; j++){
		k = i*3 + j;
		genoDat->geno[indxStrata][j] = genoRow[k];
      }
      indxStrata ++;
    }
  }

}

/****************************************************************************/

int **imatrix(int nrow, int ncol){

  int i;
  int **m;

  /* allocate pointers to rows */
  m= (int **) malloc(nrow * sizeof(int*));       
  if (!m){
    printf("Error: Allocation failure 1 in imatrix()\n");
    exit(1);
  }

  /* allocate rows and set pointers to them */
  for(i=0; i<nrow; i++){
    m[i] = (int *) malloc(ncol * sizeof(int));
    if(!m[i]){
      printf("Error: Allocation failure 2 in imatrix()\n");
      exit(1);
    }
  }


  return m;
}


/***********************************************************************************/


void printGenoDat(GENODAT *genoDat){

  int i,j;
  

  printf("\nLocus = %s, nstrata = %d:\n", genoDat->genoLabel, genoDat->nstrata);
  printf("pvalHomog = %f\n", genoDat->pvalHomog);
  printf("pHWECombo = %f\n", genoDat->phweCombined);
  printf("pHWEStrat:\n");
  for(i=0; i<genoDat->maxstrata; i++){
    printf("  p[%d] = %f\n", i+1, genoDat->phweStratum[i]);
  }


  /*
  for(i= 0; i<genoDat->nstrata; i++){
    printf("p = %f, geno:",genoDat->phweStratum[i]);
    for(j=0; j<3; j++){
      printf("%d ", genoDat->geno[i][j]);
    }
    printf("\n");
  }
  printf("\n");

*/

}


/***********************************************************************************/

void  hwdHomogExact(GENODAT * genoDat){


  int nstrata = genoDat->nstrata;
  int **geno = genoDat->geno;
  int *countA = genoDat->countA;
  int* countB = genoDat->countB;
  int *nSubj = genoDat->nSubj;
  int *hetLow = genoDat->hetLow;
  int *hetHi = genoDat->hetHi;
  int *h = genoDat->h;
  int *hLow = genoDat->hLow;
  int * hHi = genoDat->hHi;
  int* rcsumHetHi = genoDat->rcsumHetHi;
  double **lnprob = genoDat->lnprob;
 
    
  int j;
  int k;
  int m;
  int start;
  int a, het;
  int acount;
  int maxHet;

  int sumHet;
  int csum;
  int hsum;
  double lnfactObs;
  double lnfactH;

  double pvalNum = 0.0;
  double pvalDenom = 0.0;

   /* Test HWD Homogeneity Across Strata */
  if(genoDat->nstrata < 2){
	genoDat->pvalHomog = 1.0;
	return;
  }

  sumHet = 0;
  pvalNum = 0.0;
 


  for(j=0; j< genoDat->nstrata; j++){

    countA[j] = 2*geno[j][0] + geno[j][1];
    countB[j] = 2*geno[j][2] + geno[j][1];

    nSubj[j] = geno[j][0] + geno[j][1] + geno[j][2];
    sumHet += geno[j][1];
  }



  for(j=0; j< genoDat->nstrata; j++){
    hetHi[j] = ( countB[j] < countA[j] ) ? countB[j] : countA[j];
    hetLow[j] = (hetHi[j] % 2) == 0 ? 0 : 1;    
  }


  for(j=0; j< nstrata; j++){
    
    het = hetLow[j];

    /* 
       lnpprob is actually the numerator of a probability. We cannot
       compute the probability until all heterozygote configurations
       are accounted for, to divide each numerator by the sum of 
       all numerator terms.

       The next line determines lnprob for the base case, the lowest
       possible number of heterozygotes per stratum. Since this number
       is either 0 or 1, and the factorial of this is 1, we do not need
       to use het to compute lnprob for this base case.
    */


    lnprob[j][0] = lnfact(nSubj[j]) - lnfact( (countA[j] - het)/2 ) 
                 - lnfact( (countB[j] - het)/2 );

   
    for(k=1; k <= (hetHi[j]/2); k++){

      lnprob[j][k] = lnprob[j][k-1] 
	           + ((countA[j] == het) ? 0.0 : log( (double) ((countA[j] - het)/2)) )
	           + ((countB[j] == het) ? 0.0 : log( (double) ((countB[j] - het)/2)) ) 
	           - log( (double) (het + 2)) - log( (double) (het + 1));

      if(fabs(lnprob[j][k]) < TOL){
	lnprob[j][k] = 0.0;
      }


      het += 2;

    }
  }


  lnfactObs = 0.0;
  for(j=0; j<nstrata; j++){
    k = geno[j][1]/2;
    lnfactObs += lnprob[j][k];
  }

  revcumsum(hetHi, nstrata, rcsumHetHi);
  
  hLow[0] = IMAX(hetLow[0], (sumHet - rcsumHetHi[1]) );
  hHi[0]  = IMIN(sumHet, hetHi[0]);


 
  h[0] = hLow[0];
  for(j=1; j < (nstrata - 1); j++){
    csum = 0;
    for(k = 0; k < j; k++){
      csum += h[k];
    }
    
    hLow[j] = IMAX(hetLow[j], (sumHet - csum - rcsumHetHi[j+1]));
    hHi[j] = IMIN((sumHet - csum), hetHi[j]);
    h[j] = hLow[j];
  }

  hsum = 0;
  for(j=0; j< (nstrata-1); j++){
    hsum += h[j];
  }
  h[nstrata - 1] = sumHet - hsum;

  lnfactH = 0.0;
  for(j=0; j<nstrata; j++){
    lnfactH += lnprob[j][h[j]/2];
  }


  pvalDenom = exp(lnfactH);
  if( (lnfactH <= lnfactObs) || (fabs(lnfactH - lnfactObs) < TOL) ){
    pvalNum += exp(lnfactH);
  }



  while( (h[0] <= hHi[0]) && (h[nstrata-1] >= hetLow[nstrata-1])){


    if( (h[nstrata - 2] + 2) <= hHi[nstrata - 2])
      {
	h[nstrata - 2] = h[nstrata - 2] + 2;
	h[nstrata - 1] = h[nstrata - 1] - 2;
	
	lnfactH = 0.0;
	for(j=0; j<nstrata; j++){
	  lnfactH += lnprob[j][h[j]/2];
	}


	pvalDenom += exp(lnfactH);
	if( (lnfactH <= lnfactObs) || (fabs(lnfactH - lnfactObs) < TOL) ){
	  pvalNum += exp(lnfactH);
	}

      }
    else 
      {
	
	start = findIndex((nstrata - 2), h, hHi);
	if(start == (-1) ){
	  break;
	}
	
	h[start] = h[start] + 2;
	for(j=(start+1); j < (nstrata - 1); j++){
	  csum = 0;
	  for(k = 0; k < j; k++){
	    csum += h[k];
	  }
	  
	  hLow[j] = IMAX(hetLow[j], (sumHet - csum - rcsumHetHi[j+1]));
	  hHi[j] = IMIN((sumHet - csum), hetHi[j]);
	  h[j] = hLow[j];
	}
	
	hsum = 0;
	for(j=0; j< (nstrata-1); j++){
	  hsum += h[j];
	}
	h[nstrata - 1] = sumHet - hsum;



	lnfactH = 0.0;
	for(j=0; j<nstrata; j++){
	  lnfactH += lnprob[j][h[j]/2];
	}
	

	pvalDenom += exp(lnfactH);
	if( (lnfactH <= lnfactObs) || (fabs(lnfactH - lnfactObs) < TOL) ){
	  pvalNum += exp(lnfactH);
	}
	

      }
  } /* end while loop */



  genoDat->pvalHomog = exp( log(pvalNum) - log(pvalDenom));


}

/***********************************************************************************/

double lnfact(int n){
 
  double lnfact;

  if( (n == 0) || (n == 1)){
    return 0.0;
  }

  lnfact = gammln( (double)(n+1) );

  return lnfact;
}

/***********************************************************************************/

void revcumsum(int *x, int n, int *rcsum){
  int i;
  rcsum[n-1] = x[n-1];
  for(i=(n-2); i>=0; i--){
    rcsum[i] = x[i] + rcsum[i+1];
  }
}

/***********************************************************************************/

int findIndex(int i, int *h, int *hHi){
  if( (h[i] + 2) <= hHi[i] ){
    return(i);
  }

  if(i==0){
    return(-1);
  }

  return(findIndex( (i-1), h, hHi) );
}


/***********************************************************************************/
  /* Test HWE Simultaneously Across Strata */

void hweStratExact(GENODAT * genoDat){

  int nstrata = genoDat->nstrata;
  int **geno = genoDat->geno;
  int *countA = genoDat->countA;
  int * countB = genoDat->countB;
  int *nSubj = genoDat->nSubj;
  int *hetLow = genoDat->hetLow;
  int *hetHi = genoDat->hetHi;
  int *index = genoDat->index;
  int *indexUpper = genoDat->indexUpper;
  double **lnprob = genoDat->lnprob;
  double *lnpObsVec = genoDat->lnpObsVec;


  double c;
  int i;
  int j;
  int k;
  int h;
  int hobs;
  int countb;
  int stopIter;
  int countZero;

  double lnpObs;
  double sumlnprob;
  double pval;

  /* initialize pvals to equal 1.0 */
  genoDat->phweCombined = 0.0;
  for(i=0; i< genoDat->maxstrata; i++){
    genoDat->phweStratum[i] = 1.0;
  }


  for(i=0; i< nstrata; i++){

    countA[i] = 2*geno[i][0] + geno[i][1];
    countB[i] = 2*geno[i][2] + geno[i][1];

    nSubj[i] = geno[i][0] + geno[i][1] + geno[i][2];

    hetHi[i] = ( countB[i] < countA[i] ) ? countB[i] : countA[i];
    hetLow[i] = (hetHi[i] % 2) == 0 ? 0 : 1;

  }


  lnpObs  = 0.0;

  for(i=0; i<nstrata; i++){

    c = lnfact(nSubj[i]) + lnfact(countA[i]) + lnfact(countB[i]) - lnfact(2*nSubj[i]); 

    /* ln prob of observed no. heterozygotes */
    hobs = geno[i][1];

    lnpObsVec[i] =  c + hobs*log(2.) - lnfact( (countA[i] - hobs)/2 ) 
                 - lnfact(hobs) - lnfact( (nSubj[i] - (countA[i] + hobs)/2));

    lnpObs += lnpObsVec[i];

    /* each row (stratum) of lnprob is the ln prob of possible heterozygotes */
    /* base case */
    h = hetLow[i];
    lnprob[i][0] = c + h*log(2.0) - lnfact(h) 
                 - lnfact((countA[i]-h)/2) - lnfact(nSubj[i]-(countA[i]+h)/2) ;

    /* recursive cases */
    for(j=1; j<=hetHi[i]/2; j++){
      lnprob[i][j] = lnprob[i][j-1]
	           + log(4.0) + log((countA[i]-h)/2) + log(nSubj[i] - (countA[i] + h)/2)
                   - log(h + 2) - log(h + 1);
     
      h += 2;
    }

  }

 
  for(i=0; i<nstrata; i++){
    qsort(lnprob[i], (hetHi[i]/2 + 1), sizeof(double), sort_lnprob_increasing);
  }



  /* pval for hwe within each stratum */


  for(i=0; i<nstrata; i++){
    genoDat->phweStratum[i] = 0.0;
  }

  for(i=0; i<nstrata; i++){
     for(j=0; j<=hetHi[i]/2; j++){
       if( ( lnprob[i][j] - lnpObsVec[i] ) > TOL){
	 break;
       }

	 genoDat->phweStratum[i] += exp(lnprob[i][j]);

     }
  }


  for(i=0; i<nstrata; i++){
    index[i] = 0;
    indexUpper[i] = hetHi[i]/2;
  }

  pval = 0.0;
  genoDat->stopPval = 0;
  stopIter = 0;
  

  while(stopIter == 0){

    sumlnprob = 0.0;
    for(i=0; i<nstrata; i++){
      sumlnprob += lnprob[i][index[i]];
    }


    if( (sumlnprob <= lnpObs) || ( fabs(sumlnprob - lnpObs) < TOL) ){
      pval += exp(sumlnprob);
    }


    if(pval > genoDat->pStop){
      genoDat->stopPval = 1;
      stopIter = 1;
      break;
    }

    /*
      If all indices but first are equal to 0, and 
      sumlnprob > lnpObs, then incrementing any index
      will only lead to larger values of sumlnprob, so 
      stop
    */

    countZero = 0;
    for(i=1; i<nstrata; i++){
      countZero += index[i];
    }
    if( (countZero == 0) && ( (sumlnprob- lnpObs) > TOL ) ){
      stopIter = 1;
      break;
    }

    /* 
       if sumlnprob > lnpObs, then we can "jump" over the the index values that would cause
       increasing values of sumlnprob, and working backwards from last index, find the first
       index !=0, set it = 0,and increment the next index (here, next, when working
       backwords, means index[i-1]). If the next index
       exceeds its upper bound, set it = 0, and move onto the next index
       to increase. 
    */

    if( (sumlnprob - lnpObs) > TOL )
    {

      jumpIndex(index, indexUpper, nstrata, &stopIter);

      if(stopIter==1){
        break;
      }
            
    } 
    else 
      {
	/*
	  Since sumlnprob <= lnpObs, increment the appropriate index by 1, keeping
	   indices within bounds 
	*/

	incrementIndex(nstrata - 1,index, indexUpper, nstrata,  &stopIter);

	if(stopIter==1){
	    break;
	  }      
      }

  }


  genoDat->phweCombined = pval;


 
}
/***********************************************************************************/

void jumpIndex(int *index, int *indexUpper, int nstrata, int *stopIter){

  int j;

  for(j = (nstrata-1); j >= 1; j--){
    *stopIter = 1;
    if( ( index[j] > 0) | (index[j] > indexUpper[j]) ){
      *stopIter = 0;
      index[j] = 0;
      index[j-1] =  index[j-1] + 1;
      if(index[j-1] <= indexUpper[j-1]){
        return;
      }
    }
  }

  if(index[0] > indexUpper[0]){
    *stopIter = 1;
  }
  
  return;
}

/***********************************************************************************/

void incrementIndex(int j,int *index, int *indexUpper, int nstrata, int *stopIter){

  int allEqual;
  int i;

  if(*stopIter == 1){
    return;
  }


  allEqual = 1;
  for(i=0; i < nstrata; i++){
    if(index[i] != indexUpper[i])
      {
	allEqual = 0;
	break;
      }
  }

  if(allEqual == 1){
    *stopIter = 1;
    return;
  }


  index[j] = index[j] + 1;

  if( (j > 0) & (index[j] > indexUpper[j]) ){
    for(i = j; i >= 1; i--){
      if(index[i] > indexUpper[i]){
        index[i] = 0;
        index[i-1] = index[i-1] + 1;
        if(index[i-1] <= indexUpper[i-1]){
          *stopIter = 0;
          return;
        }
      }
    }
  }
  
  return;
}

/***********************************************************************************/

int sort_lnprob_increasing(const void *ptr_First, const void *ptr_Second){
 
  /* comparision function for qsort, sort on double array */

  double first =  * (double *)(ptr_First);
  double second = * (double *)(ptr_Second);

  /* for increasing sort */

  if( first  < second ) return -1;

  if( first  > second ) return  1;


  return 0;
}


/***********************************************************************************/

void hweStratHaldane(GENODAT * genoDat){

  int nstrata = genoDat->nstrata;
  int **geno = genoDat->geno;

  int i;
  int j;
  int countA;
  int countB;
  int nSubj;
  double statNum = 0.0;
  double statDenom = 0.0;
  double stat;
  double h;
  double v;

  for(i=0; i<nstrata; i++){

    h = 4.0 * geno[i][0]*geno[i][2] - geno[i][1]*(geno[i][1] - 1);
    countA = 2*geno[i][0] + geno[i][1];
    countB = 2*geno[i][2] + geno[i][1];
    nSubj = geno[i][0] + geno[i][1] + geno[i][2];

    if(nSubj > 0){
      v = 2. * countA*(countA - 1) * countB*(countB -1);
      v = v / ( (double)(2*nSubj - 3) );
      statNum   += h;
      statDenom += v;
    }

  }


  if(statDenom > TOL)
    {
      stat = (statNum * statNum) / statDenom;
      
      genoDat->pvalHaldane = 1.0 - pchisq(stat, 1);

      genoDat->zstatHaldane = statNum / sqrt(statDenom);
    }
  else{
      
    genoDat->pvalHaldane = 1.0;

    genoDat->zstatHaldane = 0.0;
  }

  return;
}

/***********************************************************************************/

void hweStratOlson(GENODAT * genoDat){
  
  int nstrata = genoDat->nstrata;
  int **geno = genoDat->geno;

  int i;
  int j;
  int countA;
  int nSubj;
  double statNum = 0.0;
  double statDenom = 0.0;
  double stat;
  double h;
  double v;
  double pA;
  double p2q2;
  int ncube;
  double num;
  double denom;
  double part1;
  double part2;


  for(i=0; i<nstrata; i++){

    h = 4.0 * geno[i][0]*geno[i][2] - geno[i][1]*(geno[i][1] - 1);
	countA = 2*geno[i][0] + geno[i][1];
    nSubj = geno[i][0] + geno[i][1] + geno[i][2];

    if(nSubj > 0){
      pA = (double)(countA) / (double) (2*nSubj);
    
      ncube = nSubj * nSubj * nSubj;

      num = (double) (8*ncube);
      denom = (double)( (2*nSubj - 1)*(2*nSubj - 2)*(2*nSubj - 3) ) ;
      part1 = (num/denom) * pA*pA*(1. -pA)*(1. -pA);

      num = (double)(2*nSubj);
      denom = (double)( (2*nSubj - 2)*(2*nSubj - 3) );
      part2 = (num/denom) * pA*(1. -pA);

      p2q2 = part1 - part2;
    
      num = (double)(nSubj*(nSubj-1));
      denom = (double)(2*(2*nSubj - 1) );

      v = (num/denom) * 4.0 * p2q2;


      statNum   += h / (double)( 2*(2*nSubj - 1) );
      statDenom += v;
    }

  }


  if(statDenom > TOL)
    {
      stat = (statNum * statNum) / statDenom;
      
      genoDat->pvalOlson = 1.0 - pchisq(stat, 1);

      genoDat->zstatOlson = statNum / sqrt(statDenom);
    }
  else{
      
    genoDat->pvalOlson = 1.0;

    genoDat->zstatOlson = 0.0;
  }

  return;
}


/***********************************************************************************/

void printHeader(int maxstrata){
  int j;

  printf("       locus          homog         hweAll");
  for(j=0; j< maxstrata; j++){
    printf("          hwe%d", j+1);
  }

  /*  
      printf("     pHaldane zHaldane");
      printf("       pOlson   zOlson");
  */

  printf("\n");
}



/***********************************************************************************/

void printHweStrat(GENODAT *genoDat){

  int j;

  printf("%12s %10.8e", genoDat->genoLabel, genoDat->pvalHomog);

  if(genoDat->stopPval == 1){
    printf(" >=%10.6e",genoDat->phweCombined);
  } else {
    printf("   %10.6e",genoDat->phweCombined);
  }

  for(j=0; j<genoDat->nstrata; j++){
    printf("  %10.6e", genoDat->phweStratum[j]);
  }

  /*
    printf(" %10.6e  %7.4f", genoDat->pvalHaldane, genoDat->zstatHaldane);
    printf(" %10.6e  %7.4f", genoDat->pvalOlson,   genoDat->zstatOlson);
  */

  printf("\n");
 
}

