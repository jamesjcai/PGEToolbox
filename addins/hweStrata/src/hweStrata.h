#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXLINE 10000
#define MAXFILENAME 25
#define MAXLABELSIZE 100
#define VERSION 0.3

/* MAX_HETEROZYGOTES is the maximum number of heterozygotes across all
   strata, used to allocate scratch memor. Setting MAX_HETEROZYGOTES to
   the maximum number of subjects in a stratum will suffice 
*/


#define MAXHETEROZYGOTES 100

#define SEPCHAR " ,\t"
#define DEBUG 0
#define TOL 0.0000001

static int imaxarg1, imaxarg2;
#define IMAX(a,b) (imaxarg1=(a), imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1, iminarg2;
#define IMIN(a,b) (iminarg1=(a), iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))


typedef struct GENODAT_T
{
  int maxstrata;
  int nstrata;
  int **geno;
  char * genoLabel;

  double pvalHomog;
  double phweCombined;
  double * phweStratum;
  int    stopPval;

  double pvalHaldane;
  double zstatHaldane;

  double pvalOlson;
  double zstatOlson;

  double pStop;


  int *countA;
  int *countB;
  int *nSubj;
  int *hetLow;
  int *hetHi;
  int *h;
  int *hLow;
  int *hHi;
  int *csumh;
  int *rcsumHetHi;
  int *index;
  int *indexUpper;
  double *lnpObsVec;

  double ** lnprob;


} GENODAT;


void initGenoDat(GENODAT *genoDat, double pStop, int maxstrata);

void printHelp(char* argvFirst, char*genoFileName, double pStop);

//unsigned long countColsFile(FILE *file);
int countColsFile(FILE *file);

void getGenoDat(char *line, int *genoRow, int maxcat, int maxstrata, GENODAT *genoDat);


int **imatrix(int nrow, int ncol);

void printGenoDat(GENODAT *genoDat);

void  hwdHomogExact(GENODAT * genoDat);

double lnfact(int n);

void revcumsum(int *x, int n, int *rcsum);

int findIndex(int i, int *h, int *hHi);

double gammln(double xx);


void hweStratExact(GENODAT * genoDat);

void jumpIndex(int *index, int *indexUpper, int nstrata, int *stopIter);

void incrementIndex(int j,int *index, int *indexUpper, int nstrata, int *stopIter);

int sort_lnprob_increasing(const void *ptr_First, const void *ptr_Second);

void hweStratHaldane(GENODAT * genoDat);

void hweStratOlson(GENODAT * genoDat);

double gammp(double a, double x);
void gcf(double *gammcf, double a, double x, double *gln);
void gser(double *gamser, double a, double x, double *gln);
double pchisq(double x, int df);
void gser(double *gamser, double a, double x, double *gln);

void printHweStrat(GENODAT *genoDat);
void printHeader(int maxstrata);
