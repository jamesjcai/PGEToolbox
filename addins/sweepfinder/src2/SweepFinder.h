#ifndef __SNPS_H__
#define __SNPS_H__


#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "factorials.h"

#define MAXPAR 200
#define BIG_ENOUGH 12.0

struct datatype {
  double loc, baselike;
  int x;
  int n;
  int folded;
} *data;

extern int datasize;
extern int nmax, nmin, xmax;
extern int sweep_width;
extern int invar;

FILE *my_fopen(char *fn, char *mode);


#endif
