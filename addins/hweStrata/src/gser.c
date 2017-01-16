/* $Author: sinnwell $ */
/* $Date: 2007/03/27 23:07:49 $ */
/* $Header: /people/biostat3/sinnwell/Projects/HWEStrat/Build/RCS/gser.c,v 1.1 2007/03/27 23:07:49 sinnwell Exp $ */
/* $Locker:  $ */
/*
 * $Log: gser.c,v $
 * Revision 1.1  2007/03/27 23:07:49  sinnwell
 * Initial revision
 * * 
 */
#include <math.h>
#include <stdio.h>

#define ITMAX 100
#define EPS 3.0e-7

void gser(double *gamser, double a, double x, double *gln)
{
	double gammln(double xx);
	void nrerror(char error_text[]);
	int n;
	double sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) printf("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		printf("a too large, ITMAX too small in routine gser");
		return;
	}
}
#undef ITMAX
#undef EPS
