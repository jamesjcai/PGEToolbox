/* $Author: sinnwell $ */
/* $Date: 2007/03/27 23:07:35 $ */
/* $Header: /people/biostat3/sinnwell/Projects/HWEStrat/Build/RCS/get_optval.c,v 1.1 2007/03/27 23:07:35 sinnwell Exp $ */
/* $Locker:  $ */
/*
 * $Log: get_optval.c,v $
 * Revision 1.1  2007/03/27 23:07:35  sinnwell
 * Initial revision
 * * 
 */

/*
  Get options and their valuues from string, where the string was
  processed by cmdString.

  An option is identified by starting with a '-', and its 
  assigned value is the value that follows a matched option.
  
  The type of value is identified by vtype (see header file for
  possible types vtypes).
  
  Returned value = 0 if option not found
                   1 if option found and its value found
                  -1 if option found but its value not found

  Example usage in a main program as follows.




  int help=0;
  double pStop = 1.0;
  char *genoFileName = "geno.dat";
  int retOptVal;

  char * cmd = cmdString(argc, argv);

  if(get_optval(cmd,"h",OV_FLAG,&help) == 1) 
    {
      printHelp(argv[0], genoFileName, pStop);
      exit(1);
    }

  if( get_optval(cmd, "geno",OV_STRING, genoFileName) == -1)
    {
      printf("Error: no value for option geno\n");
      exit(1);
    }

  if( get_optval(cmd,"pstop",OV_DOUBLE, &pStop) == -1)
    {
      printf("Error: no value for option pstop\n");
      exit(1);
	}


*/



#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "get_optval.h"

int get_optval(char * argv, char *opt, int vtype, void *val) {

   int found, j, k, len, iv, negate, argnum;
   float fv;
   long  lv;
   unsigned int uiv;
   unsigned long ulv;
   double dv;
   char *strv;

   char fs[MFLEN];
   char *argi;

   found = 0;

   len = strlen(opt);
   if (!len) return 0;

   argi = argv;

   while(*argi != '\0'){


     if (*argi == '-') {
 
       if (argi[1]=='n' && argi[2]=='o') 
	 {
	   argi += 3;
	   negate = 1;
	 }
       else 
	 {
	   argi ++;
	   negate = 0;
	 }


       found = 1;

       for (j=0; j<len; j++)
	 {
	   found = found && (*argi == opt[j]);
	   argi ++;
	 }

     }
    
     if (found==1) 
       {
	 if (negate && vtype)
	   {
	     return -1; 
	   }

	 k = 0;

	 /* 
	    to allow a negative value, need to allow for 
	    the first char of a value to be '-'
	 */
	 if( (vtype >=1) && (vtype <=4) && (*argi == '-') )
	   {
	     fs[k] = *argi;
	     argi ++;
	     k ++;
	   }

	 while( (*argi != '-') && (*argi != '\0') )
	   {

	     /* skip over equal sign */

	     if(*argi == '='){
	       argi ++;
	     }

	     fs[k] = *argi;	 
	     argi ++;
	     k ++;
	   }

	 fs[k] = '\0';

	 if( (vtype > 0) && (strlen(fs)==0) )
	   {
	     return -1;
	   }

	 switch (vtype)
	   {

	   case 0:
	     if (strlen(fs)>1) return -1;
	     if (strlen(fs)==0) {
	       *(int*)val = negate ? 0 : 1;
	       return 1;
	     }

	     if (fs[0] ==  '+') {
	       *(int*)val = negate ? 0 : 1;
	       return 1;
	     }

	     if (fs[0] ==  '-') {
	       *(int*)val = negate ? 1 : 0;
	       return 1;
	     }

	     return -1;

	   case 1:
	     if (sscanf(fs, "%d", &iv)) 
	       {
		 *(int*)val = iv;
		 return 1;
	       } 
	     else 
	       {
		 return -1;
	       }

	   case 2:
	     if (sscanf(fs, "%d", &lv)) 
	       {
		 *(long*)val = lv;
		 return 1;
	       } 
	     else 
	       {
		 return -1;
	       }

	   case 3:
	     if (sscanf(fs, "%f", &fv)) 
	       {
		 *(float*)val = fv;
		 return 1;
	       } 
	     else 
	       {
		 return -1;
	       }

	   case 4:
	     if (sscanf(fs, "%lf", &dv)) 
	       {
		 *(double*)val = dv;
		 return 1;
	       } 
	     else 
	       {
		 return -1;
	       }

	   case 5:
	     if (sscanf(fs, "%d", &uiv)) 
	       {
		 *(unsigned int*)val = uiv;
		 return 1;
	       } 
	     else 
	       {
		 return -1;
	       }

	   case 6:
	     if (sscanf(fs, "%d", &ulv)) 
	       {
		 *(unsigned long*)val = ulv;
		 return 1;
	       } 
	     else 
	       {
		 return -1;
	       }

	   case 7:

	     strv = (char *) malloc( (1 + strlen(fs)) * sizeof(char) );

	     strcpy( strv, fs);

	     *(char **) val = strv;

	     return 1;

	   }
       }
     else 
       {
	 argi ++;
       }
   }


   return found; /* Opt not found */

}

/******************************************************************************/

char * cmdString(int argc, char ** argv){

  int i, j, k, len;
  int totalLength = 0;
  char *cmdString;

  for(i=1; i < argc; i++){
	totalLength += strlen(argv[i]);
  }

  cmdString = (char *) malloc((totalLength + 1) * sizeof(char) );

  k = 0;

  for(i=1; i < argc; i++){
    len = strlen(argv[i]);
	for(j=0; j<len; j++){
      if(argv[i][j] == '>'){
		k++;
		break;
	  }

      cmdString[k] = argv[i][j];
      k++;
    }
  }

  cmdString[k] = '\0';

  return(cmdString);
}
