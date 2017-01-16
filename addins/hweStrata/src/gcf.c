/* $Author: sinnwell $ */
/* $Date: 2007/03/27 23:07:26 $ */
/* $Header: /people/biostat3/sinnwell/Projects/HWEStrat/Build/RCS/gcf.c,v 1.1 2007/03/27 23:07:26 sinnwell Exp $ */
/* $Locker:  $ */
/*
 * $Log: gcf.c,v $
 * Revision 1.1  2007/03/27 23:07:26  sinnwell
 * Initial revision
 * * 
 */
#include <math.h>
#include <stdio.h>

#define ITMAX 500
#define EPS 3.0e-7
#define FPMIN 1.0e-30

void gcf(double *gammcf, double a, double x, double *gln)
{
	double gammln(double xx);
	void nrerror(char error_text[]);
	int i;
	double an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) printf("a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN
