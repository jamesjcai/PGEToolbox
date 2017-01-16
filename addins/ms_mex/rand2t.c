/*  Link in this file for random number generation with rand()
     and seeding from the clock  */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double
ran1()
{
	int rand();	
	return( rand()/(RAND_MAX+1.0)  );
}


	void seedit( const char *flag)
{	
	unsigned int seed2 ;

  if( flag[0] == 's' ) {
	/* srand( seed2 = time(NULL)) ; */
	  srand( seed2 = time(0)) ;
    /* mexPrintf("\n%d\n", seed2 ); */
}
}
