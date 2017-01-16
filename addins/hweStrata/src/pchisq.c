/* $Author: sinnwell $ */
/* $Date: 2007/03/27 23:07:17 $ */
/* $Header: /people/biostat3/sinnwell/Projects/HWEStrat/Build/RCS/pchisq.c,v 1.1 2007/03/27 23:07:17 sinnwell Exp $ */
/* $Locker:  $ */
/*
 * $Log: pchisq.c,v $
 * Revision 1.1  2007/03/27 23:07:17  sinnwell
 * Initial revision
 * * 
 */

double gammp(double a, double x);

double pchisq(double x, int df){

  /* return cdf P(X <=x) for chi-square with df */
  
  double p;

  p = gammp( ((double) df )/2.0 , x/2.0);

  return p;

}
