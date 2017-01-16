/*  Link in this file for random number generation with rand() */

#include <stdio.h>
#include <stdlib.h>


	double
ran1()
{
	int rand();
    double x;
    x=rand()/(RAND_MAX+1.0);
    //if(DEBUG){mexPrintf("%f ",x); mexPrintf("\n");}
	return( x  );
}

	void seedit( const char *flag)
{
	FILE *fopen(), *pfseed;
	unsigned int seed2 ;
  if( flag[0] == 's' ) {
    pfseed = fopen("seedms","r");
        if( pfseed == NULL ) {
           seed2 = 59243; }
        else {
          fscanf(pfseed," %d",&seed2);
          fclose( pfseed);
          }
          srand( seed2) ;

        mexPrintf("\n%d\n", seed2 );    
	}
   else {
    	pfseed = fopen("seedms","w");
        fprintf(pfseed,"%d \n",rand());  
        fclose( pfseed);
	}
}

	int
commandlineseed( char **seeds)
{
	unsigned int seed2 ;
    void srand(unsigned int seed);
	
	seed2 = atoi( seeds[0] );

	printf("\n%d\n", seed2 );    

	srand(seed2) ; 
	return(1);
}
