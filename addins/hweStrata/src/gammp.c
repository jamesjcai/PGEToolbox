/* $Author: sinnwell $ */
/* $Date: 2007/03/27 23:06:57 $ */
/* $Header: /people/biostat3/sinnwell/Projects/HWEStrat/Build/RCS/gammp.c,v 1.1 2007/03/27 23:06:57 sinnwell Exp $ */
/* $Locker:  $ */
/*
 * $Log: gammp.c,v $
 * Revision 1.1  2007/03/27 23:06:57  sinnwell
 * Initial revision
 * * 
 */
#include <stdio.h>

double gammp(double a, double x)
{
	void gcf(double *gammcf, double a, double x, double *gln);
	void gser(double *gamser, double a, double x, double *gln);
	void nrerror(char error_text[]);
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) {
	  printf("Invalid arguments in routine gammp");
	}

	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}
