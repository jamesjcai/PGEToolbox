#include <stdarg.h>
#include <math.h>
#include "SweepFinder.h"
#include "my_rand.h"
#include "freq.h"
#include "sort.h"
//#include "matrix.h"

//stores the probabilities over a grid of values of alpha*d
//indexed by prob[n][x][2][gridsize]
//gridsize is hard-coded, not same as argument to program
//prob[n][x][0][g] holds probabilities, 
//prob[n][x][1][g] holds vector that spline function needs
double ****prob=NULL;
double *ads=NULL; //holds vector of alpha*d values used to get prob

int invar=0;
int datasize=0, sweep_width, nmax, nmin, xmax;
double lowbound[MAXPAR], upbound[MAXPAR];
struct datatype *data=NULL;

void getlims() {
  int i;
  invar=0;
  if (datasize==0) {
    nmin=nmax=0;
    return;
  }
  nmin = nmax = data[0].n;
  for (i=0; i<datasize; i++) {
    if (data[i].n>nmax) nmax=data[i].n;
    if (data[i].n<nmin) nmin=data[i].n;
    if (invar==0 && data[i].x==data[i].n) invar=1;
    if (invar!=2 && data[i].x==0 && data[i].n > 0) invar=2;
  }
  if (invar) xmax=nmax+1;
  else xmax=nmax;
}

FILE *my_fopen(char *fn, char *mode) {
  FILE *rv = fopen(fn, mode);
  if (rv==NULL) {
    fprintf(stderr, "error opening %s mode \"%s\"\n", fn, mode);
    exit(-1);
  }
  return rv;
}


void get_nlim() {
  int i;
  if (datasize==0) {
    nmax=nmin=0;
    return;
  }
  nmax=nmin=data[0].n;
  for (i=1; i<datasize; i++) {
    if (data[i].n > nmax) nmax=data[i].n;
    if (data[i].n < nmin) nmin=data[i].n;
  }
}


void get_lims() {
  fflush(stdout);
  get_nlim();
  if (invar) xmax = nmax+1;
  else xmax = nmax;
}

#define NCOLUMN 4
#define LOC 0
#define X 1
#define N 2
#define FOLDED 3

//all column headers except folded are required.
//if folded is not present, all SNPs will be assumed unfolded
char colname[NCOLUMN][10]={"position", "x", "n", "folded"};

//get rid of some data so we can test missing data
void losedata(double per) {
  struct datatype *newdata = malloc(datasize*sizeof(struct datatype));
  int newdatasize=0;
  int i,j , x, n;
  for (i=0; i<datasize; i++) {
    x = data[i].x;
    n = data[i].n;
    for (j=0; j<data[i].n; j++) {
      if (uniform()  < per) {
	if (j < data[i].x) x--;
	n--;
      }
    }
    if (x!=0 && x!=n) {
      newdata[newdatasize].x = x;
      newdata[newdatasize].n = n;
      newdata[newdatasize].folded=data[i].folded;
      newdata[newdatasize].loc = data[i].loc;
      newdatasize++;
    }
  }
  if (datasize > 0) free(data);
  datasize = newdatasize;
  data = newdata;
  data = realloc(newdata, datasize*sizeof(struct datatype));
  printf("done losedata datasize=%i\n", datasize);
}

void readms_error(char *infn) {
  printf("error reading data from %s\n", infn);
  exit(-1);
}


int readsnps_ms(char *infn) {
  static FILE *infile=NULL;
  static char *curr_infile=NULL;
  static int numchr=0, rep=0;
  char str[1000], c;
  int i, j;

  if (infn!=curr_infile) {
    if (infile!=NULL) 
      fclose(infile);
    printf("opening %s\n", infn);
    rep=0;
    infile=fopen(infn, "r");
    curr_infile=infn;
    if (EOF==fscanf(infile, "%s %i", str, &numchr)) readms_error(infn);
    if (strstr(str, "ms")==NULL || numchr <=0) readms_error(infn);
    printf("numchr=%i\n", numchr);
  }
  while (1) {
    if (EOF==fscanf(infile, "%s", str)) {
      fclose(infile);
      infile=NULL;
      curr_infile=NULL;
      return 0;
    }
    if (strcmp(str, "segsites:")==0) break;
  }
  if (datasize > 0) free(data);
  datasize=0;
  if(EOF==fscanf(infile, "%i",&datasize)) readms_error(infn);
  printf("reading %s rep %i\n", infn, ++rep);
  printf("datasize=%i\n", datasize);
  if (datasize<=0) {
    printf("readsnps_ms: segsites must be greater than 0!\n");
    readms_error(infn);
  }
  data = malloc(datasize*sizeof(struct datatype));
  if(EOF==fscanf(infile, "%s", str)) readms_error(infn);
  if(strcmp(str, "positions:")!=0) readms_error(infn);
  for (i=0; i<datasize; i++) {
    if (EOF==fscanf(infile, "%lf", &data[i].loc))
      readms_error(infn);
    data[i].folded=0;
    data[i].x=0;
    data[i].n=numchr;
  }
  for (i=0; i<numchr; i++) {
    while ('\n'!=(c=fgetc(infile))) if (isspace(c)==0) readms_error(infn);
    for (j=0; j<datasize; j++) {
      c=fgetc(infile);
      if (c!='0' && c!='1') readms_error(infn);
      data[j].x += (c=='1');
    }
  }
  for (i=0; i<datasize; i++) {
    if (data[i].x==0 || data[i].x>=data[i].n) {
      printf("invalid data in msfile %s data[%i].x=%i n=%i\n", infn,
	     i, data[i].x, data[i].n);
      readms_error(infn);
    }
  }
  getlims();
  return 1;
}

  
  
      
      


//this is the format the old version of yuseob's program outputs
//it doesn't output n anywhere so you have to know what it is
void readsnps_yuseob(char *infn, int rep, int n) {
  FILE *infile = my_fopen(infn, "r");
  char str[1000], c;
  int count=0, nsites, i;
  if (datasize > 0) free(data);
  datasize=0;
  while (1) {
    while (EOF!=fscanf(infile, "%s", str))
      if (strcmp(str, "sites")==0) break;
    if (strcmp(str, "sites")!=0) break;
    while (':'!=(c=fgetc(infile))) assert(isspace(c));
    fscanf(infile, "%i", &nsites);
    if (rep < 0 || count==rep) {
      if (datasize==0) data = malloc(nsites*sizeof(struct datatype));
      else data = realloc(data, (datasize+nsites)*sizeof(struct datatype));
      for (i=0; i<nsites; i++) {
	assert(EOF!=fscanf(infile, "%lf %i", &data[datasize+i].loc,
			   &data[datasize+i].x));
	data[datasize+i].folded=0;
	data[datasize+i].n=n;
      }
      datasize += nsites;
      if (count==rep) break;
    }
    count++;
  }
  fclose(infile);
  losedata(0.1);
  getlims();
  printf("done getlims_yuseob datasize=%i nmax=%i nmin=%i xmax=%i invar=%i\n",
	 datasize, nmax, nmin, xmax, invar);
}


void printdata(char *outfn) {
  int i;
  FILE *outfile = fopen(outfn, "w");
  fprintf(outfile, "position\tx\tn\tfolded\n");
  for (i=0; i<datasize; i++)
    fprintf(outfile, "%f\t%i\t%i\t%i\n", data[i].loc,
	    data[i].x, data[i].n, data[i].folded);
  fclose(outfile);
}


void readsnps_error(char *infn) {
  printf("error reading %s\n", infn);
  exit(-1);
}


void readsnps(char *infn) {
  FILE *infile=my_fopen(infn, "r");
  int colpos[NCOLUMN], pos=0, col=0, i, j;
  char c, str[1000];
  if (datasize > 0) free(data);
  datasize=0;
  c=fgetc(infile);
  for (i=0; i<NCOLUMN; i++) colpos[i]=-1;
  while (c!='\n' && c!=EOF) {
    str[0]=c;
    pos=1;
    while ('\t'!=(c=fgetc(infile)) && c!='\n' && c!=EOF)
      str[pos++]=c;
    str[pos]='\0';
    for (i=0; i<NCOLUMN; i++)
      if (strcmp(str, colname[i])==0) {
	colpos[i]=col;
	break;
      }
    col++;
    if (c=='\n' || c==EOF) break;
    c=fgetc(infile);
  }
  if (colpos[LOC]<0 || colpos[X]<0 || colpos[N]<0) {
    fprintf(stderr, "readsnps: infile should have columns named position, x, and n (and optionally folded\n");
    exit(-1);
  }
  while (EOF!=(c=fgetc(infile)))
    if (c=='\n') datasize++;
  fclose(infile);
  infile=my_fopen(infn, "r");
  while ('\n'!=(c=fgetc(infile)) && c!=EOF);
  data = malloc(datasize*sizeof(struct datatype));
  for (i=0; i<datasize; i++) {
    for (j=0; j<col; j++) {
      if (colpos[LOC]==j) 
	{if(EOF==fscanf(infile, "%lf", &data[i].loc)) readsnps_error(infn);}
      else if (colpos[X]==j) 
	{if(EOF==fscanf(infile, "%i", &data[i].x)) readsnps_error(infn);}
      else if (colpos[N]==j) 
	{if(EOF==fscanf(infile, "%i", &data[i].n)) readsnps_error(infn);}
      else if (colpos[FOLDED]==j) 
	{if(EOF==fscanf(infile, "%i", &data[i].folded)) readsnps_error(infn);}
      else {
	while ('\t'!=(c=fgetc(infile)) && c!='\n' && c!=EOF);
	if (c=='\n' || c==EOF) if(j!=col-1) readsnps_error(infn);
      }
    }
  }
  fclose(infile);
  if (colpos[FOLDED] < 0)
    for (i=0; i<datasize; i++) 
      data[i].folded=0;
  getlims();
  printf("done readsnps datasize=%i nmax=%i nmin=%i xmax=%i invar=%i\n", datasize, nmax, nmin, xmax, invar);
}



//probability of choosing j from sample size s given frequency 
//spectrum p
double p_j_s(int j, int s, double *p) {
  int i;
  double rv=0.0;
  if (j < 0 || j > s) return 0.0;
  for (i=j; i<=nmax; i++)
    if (s-j >=0 && i<xmax)
      rv += p[i]*xchoosey(i,j)*xchoosey(nmax-i,s-j)/
	xchoosey(nmax,s);
  return rv;
}

	 
double p_kescapes(int k, double ad) {
  double pe = exp(-ad);
  return xchoosey(nmax, k)*pow((1.0-pe),k)*pow(pe,(nmax-k));
}


double get_pstar_sub(int b0, double* p, double ad) {
  int k, b;
  static double lastad=-1.0;
  static double *rv=NULL;

  if (lastad!=ad) {
    if (rv==NULL) rv = malloc((nmax+1)*sizeof(double));
    for (b=0; b<=nmax; b++) {
      rv[b] = 0.0;
      if (b < xmax)
	rv[b] = p[b]*p_kescapes(nmax, ad);
      for (k=0; k<nmax; k++) {
	rv[b] += p_kescapes(k, ad)*
	  (p_j_s(b+1-nmax+k, k+1, p)*(b+1-nmax+k)/(k+1) +
	   p_j_s(b, k+1, p)*(k+1-b)/(k+1));
      }
    }
    lastad = ad;
  }
  return rv[b0];
}

double get_pstar(int i, double* p, double ad) {
  double pr, sum;

  sum=1.0;
  if (invar==0 || invar==1)
    sum -= get_pstar_sub(0, p, ad);
  if (invar==0)
    sum -= get_pstar_sub(nmax, p, ad);
  if (sum > 1.0 || sum <= 0.0) {
    int j;
    sum = 0.0;
    for (j=1; j<nmax; j++) sum += get_pstar_sub(j, p, ad);
    if (invar==1) sum += get_pstar_sub(nmax, p, ad);
  }
  assert(sum <= 1.0 && sum > 0.0);
  pr = get_pstar_sub(i, p, ad)/sum;
  return pr;
}



/*double get_pstar(int i, double* p, double ad) {
  double pr, sum;

  sum=1.0;
  if (invar==0 || invar==1)
    sum -= get_pstar_sub(0, p, ad);
  if (invar==0) 
    sum -= get_pstar_sub(nmax, p, ad);
  

  assert(sum <= 1.0 && sum > 0.0);
  pr = get_pstar_sub(i, p, ad)/sum;
  return pr;
  }*/

double splint(double* xvec, double *yvec, double *yvec2, int n, 
	      double x)
{
  int lowpos,hipos,pos;
  double diff,b,a;
  double y;
  
  if (x < xvec[0]) return xvec[0];
  lowpos=0;
  hipos=n-1;
  while (hipos-lowpos > 1) {
    pos=(hipos+lowpos)/2;
    if (xvec[pos] > x) hipos=pos;
    else lowpos=pos;
  }
  diff=xvec[hipos]-xvec[lowpos];
  if (diff == 0.0) {fprintf(stderr, "splint error\n"); exit(-1);}
  a=(xvec[hipos]-x)/diff;
  b=(x-xvec[lowpos])/diff;
  y=a*yvec[lowpos]+b*yvec[hipos]+((a*a*a-a)*yvec2[lowpos]+(b*b*b-b)*yvec2[hipos])*(diff*diff)/6.0;
  if (y<=0.0 || y>=1.0) {
    y=(yvec[lowpos]*a+b*yvec[hipos]);
  }
  if (y <0.0 || y>1.0) {
    printf("y=%e lowpos=%i hipos=%i a=%e b=%e yvec=%e, %e\n", 
	   y, lowpos, hipos, a, b, yvec[lowpos], yvec[hipos]);
    y = (yvec[lowpos]+yvec[hipos])/2.0;
  }
  return y;
}


void spline(double *x, double *y ,int n, double yp1, double ypn,
	    double *y2)
{
  int i,k;
  double p,qn,sig,un,*u;
  
  /*	for (i=1; i<=n; i++) {
	
  printf("%i\t%e\t%e\t%e\n", i, x[i], y[i], y2[i]);
  }
  exit(-1);*/
  
  u=malloc(n*sizeof(double));
  if (yp1 > 0.99e30)
    y2[0]=u[0]=0.0;
  else {
    y2[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
    if (x[1]-x[0]==0.0) {
      printf("spline: x[1]=x[0]=%e\n", x[1]); exit(-1);}
    
  }
  for (i=1;i<=n-2;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    if (x[i+1]==x[i-1]) {
      printf("spline2: x[%i]=x[%i]=%e\n", i+1, i-1, x[i+1]);
      exit(-1);}
    p=sig*y2[i-1]+2.0;
    if (p==0.0) {
      printf("spline: p=%e\n", p); exit(-1);}
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    if (x[i+1]==x[i]) {
      printf("spline3: x[%i]=x[%i]=%e\n", i+1, i, x[i+1]); exit(-1);}
    if (x[i]==x[i-1]) {
      printf("spline4: x[%i]=x[%i]=%e\n", i, i-1, x[i]); exit(-1);}
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    printf("spline: Else!\n"); exit(-1);
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  if (qn*y2[n-2]+1.0==0.0) {
    printf("spline: qn=%e y2[%i]=%e\n", qn, n-2, y2[n-2]); exit(-1);}
  for (k=n-2;k>=0;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  free(u);
}


void calcprobs(double *p, int gridsize) {
  double pr, maxval, minval, interval;
  int i, n, x, j, k, minx, maxx, maxxi, xi;
  printf("calcprobs nmax=%i nmin=%i xmax=%i invar=%i\n",
	 nmax, nmin, xmax, invar);
  fflush(stdout);

  ads = malloc(gridsize*sizeof(double));
  //gridsize=60 scheme
  //this is a bit ad-hoc but i'm choosing values of ad that cover
  //the entire range of minval to maxval, but concentrate on the higher
  //values which are usually more relevant
  //so first 20 values are chosen on a log scale, other 40 on linear
  minval=1.0e-8;
  maxval=BIG_ENOUGH;
  interval=log(maxval/minval)/19.0;
  for (i=0; i<20; i++)
    ads[i]=exp(log(minval)+i*interval);
  minval=0.1;
  maxval=3.0;
  interval=(maxval-minval)/39.0;
  for (i=0; i<40; i++)
    ads[i+20]=minval+interval*i;
  dsort(gridsize, ads-1);
  for (i=1; i<gridsize; i++)
    if (ads[i]==ads[i-1]) {
      printf("ads[%i]=ads[%i]=%e\n\n\n", i, i-1, ads[i]);
      exit(-1);
    }

  prob = malloc((nmax+1)*sizeof(double***));
  for (n=nmin; n<=nmax; n++) {
    if (invar==0) {minx=1; maxx=n;}
    else if (invar==1) {minx=1; maxx=n+1;}
    else {minx=0; maxx=n+1;}
    prob[n] = malloc(maxx*sizeof(double**));
    for (x=minx; x<maxx; x++) {
      prob[n][x] = malloc(2*sizeof(double*));
      for (j=0; j<2; j++) {
	prob[n][x][j] = malloc(gridsize*sizeof(double));
	for (k=0; k<gridsize; k++) 
	  prob[n][x][j][k] = 0.0;
      }
    }
  }
  for (k=0; k<gridsize; k++) {
    fflush(stdout);
    for (n=nmin; n<=nmax; n++) {
      if (invar==0) {minx=1; maxx=n;}
      else if (invar==1) {minx=1; maxx=n+1;}
      else {minx=0; maxx=n+1;}
      for (x=minx; x<maxx; x++) {
	maxxi = x + nmax - n;
	for (xi=x; xi<=maxxi; xi++) {
	  pr = exp(ln_xchoosey(xi, x)+
		   ln_xchoosey(nmax-xi, n-x)-
		   ln_xchoosey(nmax, n));
	  prob[n][x][0][k] +=
	    pr*get_pstar(xi, p, ads[k]);
	}
      }
    }
  }
  
  for (n=nmin; n<=nmax; n++) {
    if (invar==0) {minx=1; maxx=n;}
    else if (invar==1) {minx=1; maxx=n+1;}
    else {minx=0; maxx=n+1;}
    for (x=minx; x<maxx; x++) {
      spline(ads, prob[n][x][0], gridsize, 1.0e30, 
	     1.0e30, prob[n][x][1]);
    }
  }
  printf("done calcprob\n");
  fflush(stdout);
}


//actually returns LR now
double ln_likelihood(double *p, double alpha, double sweep) {
  int s, n, x, gridsize=60, first;
  double like=0.0, dist, ad, pr;

  if (prob==NULL)
    calcprobs(p, gridsize);

  sweep_width=0;
  for (s=0; s<datasize; s++) {
    dist = fabs(data[s].loc-sweep);
    if (dist < 0.01) dist = 0.01;  //don't want to be right on top of sweep
    ad = alpha*dist;
    if (ad >= BIG_ENOUGH) continue;
    sweep_width++;
    n=data[s].n;
    x = data[s].x;
    pr =0;
    first=1;
  calc:
    if (ad < ads[0]) 
      pr += prob[n][x][0][0];
    else 
      pr += splint(ads, prob[n][x][0],
		   prob[n][x][1], 
		   gridsize, ad);
    
    if (pr >= 0.0 && pr < 1.0);
    else {
      printf("1: val=%e\n", pr);
      printf("%i %i %e\n", nmax-n, x, ad);
      printf("%e %e\t%e %e\n", ads[gridsize-2], 
	     prob[n][x][0][gridsize-2],
	     ads[gridsize-1],
	     prob[n][x][0][gridsize-1]);
      //      debug_splint=1;
      splint(ads, prob[n][x][0],
	     prob[n][x][1],
	     gridsize, ad);
      
    }    
    if (pr < 0.0 || pr >= 1.00000001) {
      printf("val=%e\n", pr); exit(-1);
    }
    assert(pr >=0.0 && pr<1.00000001);
    
    if (first && data[s].folded) {
      x = n-x;
      first=0;
      if (x!=n-x)
	goto calc;
    }
    like += log(pr) - data[s].baselike;
  }
  return like;
}

//returns maximum likelihood for sweep at given location (sweep) 
//and also sets MLE for alpha
double findmax_alpha(double *alpha, double *pi, double sweep) {
  int startgridsize=100, gridsize, max=-1, i, sw[100], count=0;
  double *vals, *likes, tol=1.0e-6, minalpha, maxalpha, mind=-1, maxd=-1, 
    interval, like, totd=0;
  gridsize=startgridsize;
  vals = malloc(gridsize*sizeof(double));
  likes = malloc(gridsize*sizeof(double));
  for (i=0; i<datasize; i++) {
    if (i==0 || fabs(data[i].loc-sweep) < mind)
      mind = fabs(data[i].loc-sweep);
    if (i==0 || fabs(data[i].loc-sweep) > maxd)
      maxd = fabs(data[i].loc-sweep);
    totd += fabs(data[i].loc-sweep);
  }
  if (mind<0.01) mind=0.01;
  maxalpha=(BIG_ENOUGH+1)/mind;
  minalpha = BIG_ENOUGH/(totd/datasize);
  //  printf("minalpha=%e maxalpha=%e\n", minalpha, maxalpha);

  while ((maxalpha-minalpha)/((maxalpha+minalpha)/2.0) > tol) {
    interval = log(maxalpha/minalpha)/(gridsize-1);
    for (i=0; i<gridsize; i++) {
      vals[i] = exp(log(minalpha)+i*interval);
      if (i!=0 && sw[i-1]==0)
	likes[i] = likes[i-1];
      else likes[i] = ln_likelihood(pi, vals[i], sweep);
      //      printf("%i %f %e like=%f sweep_width=%i\n", i, sweep, vals[i], likes[i], sweep_width);
      sw[i] = sweep_width;
      if (i==0 || likes[i] > likes[max])
	max = i;
    }
    if (max==0) 
      minalpha = exp(log(minalpha)-gridsize*interval);
    else minalpha = vals[max-1];
    if (max==gridsize-1) maxalpha += (vals[max-1]-vals[max-2]);
    else maxalpha = vals[max+1];
    gridsize = 5;
    count++;
    /*    if (count > 10000) {
      printf("%i %e %e\n", count, minalpha, maxalpha);
      fflush(stdout);
      }*/
  }
  like = likes[max];
  sweep_width = sw[max];
  *alpha = vals[max];
  free(vals);
  free(likes);
  return like;
}


double findsweeps(char *outfn, double *p, int gridsize, int noisy, int msrep) {
  double alpha, lr, maxlr=0.0, minlike=0.0, smax=0, smin=1e100, sweep;
  int i, rep;
  FILE* outfile;
  double maxalpha=-1.0, maxsweep=-1.0;

  for (i=0; i<datasize; i++) 
    minlike += (data[i].baselike = 
		likelihood_freq_onesnp(data[i].x, data[i].n,
				       data[i].folded, 1, p));
  for (i=0; i<datasize; i++) {
    if (data[i].loc < smin || i==0) smin = data[i].loc;
    if (data[i].loc > smax || i==0) smax = data[i].loc;
  }
  printf("findsweeps smin=%e smax=%e gridsize=%i minlike=%f\n", 
	 smin, smax, gridsize, minlike);

  if (noisy) {
    outfile=my_fopen(outfn, "w");
    fprintf(outfile, "location\tLR\talpha\n");
    fclose(outfile);
  }
  if (noisy==0 && msrep==1) {
    outfile = my_fopen(outfn, "w");
    fprintf(outfile, "rep\tLR\tlocation\talpha\n");
    fclose(outfile);
  }

  for (rep=0; rep<gridsize; rep++) {
    if (gridsize==1) sweep=(smin+smax)/2.0;  //gridsize=1 is not recommended
    else sweep = smin + (smax-smin)*rep/(gridsize-1);
    lr = findmax_alpha(&alpha, p, sweep);
    if (rep==0 || lr > maxlr) {
      maxlr = lr;
      maxalpha = alpha;
      maxsweep = sweep;
    }
    if (noisy) {
      outfile = my_fopen(outfn, "a");
      fprintf(outfile, "%f\t%f\t%e\n", sweep, lr, alpha);
      fclose(outfile);
      printf("pos %f\tLR=%f\talpha=%e\n", sweep, lr, alpha);
    }
  }
  if (noisy==0) {
    outfile=my_fopen(outfn, "a");
    fprintf(outfile, "%i\t%f\t%f\t%e\n", msrep, maxlr, maxsweep, maxalpha);
    fclose(outfile);
  }
  printf("maxsweep LR=%f loc=%f alpha=%e", maxlr, maxsweep, maxalpha);
  if (msrep!=-1) printf("\tmsrep=%i\n", msrep);
  else fputc('\n', stdout);
  return maxlr;
}


void usage() {
  printf("usage:\n");
  printf("\tto get frequency spectrum: ./SweepFinder -f infilename outfilename \n");
  printf("\tto find sweeps: ./SweepFinder -s GRIDSIZE SNPfilename outfilename\n");
  printf("\tto analyze ms output with optimized frequency spectrum for each rep: ./SweepFinder -m GRIDSIZE msfilename outfilename\n");
  printf("\tto analyze ms output with same frequency spectrum for each rep: ./Sweepfinder -M GRIDSIZE msfilename freqfilename outfilename\n");
  printf("\tto find sweeps using a pre-computed frequency spectrum: ./SweepFinder -l GRIDSIZE infilename freqfilename outfilename\n");
  exit(-1);
}


int main(int argc, char *argv[]) {
  double *p;
  int gridsize, rep;
  char snpfn[1000], outfn[1000], freqfn[1000];
  if (argc < 3) usage();
  if (strlen(argv[1])!=2 || argv[1][0]!='-')
    usage();
  if (argv[1][1]=='f') {
    if (argc!=4) usage();
    sprintf(snpfn, argv[2]);
    sprintf(outfn, argv[3]);
    readsnps(snpfn);
    //readsnps_yuseob(snpfn, -1, 50);
    getfreq(outfn);
  }
  else if (argv[1][1]=='l') {
    if (argc!=6) usage();
    gridsize = atoi(argv[2]);
    if (gridsize <=0) {
      printf("gridsize should be > 0!\n"); usage();
    }
    sprintf(snpfn, argv[3]);
    sprintf(freqfn, argv[4]);
    sprintf(outfn, argv[5]);
    readsnps(snpfn);
    //    readsnps_yuseob(snpfn, 2, 50);
    //    printdata("testsnps.txt");
    //    exit(-1);
    p = loadfreq(freqfn);
    findsweeps(outfn, p, gridsize, 1, -1);
    free(p);
  }
  else if (argv[1][1]=='s') {
    if (argc!=5) usage();
    gridsize = atoi(argv[2]);
    if (gridsize<=0) {
      printf("gridsize should be > 0\n"); usage();
    }
    sprintf(snpfn, argv[3]);
    sprintf(freqfn, "tempfreq.txt");
    sprintf(outfn, argv[4]);
    readsnps(snpfn);
    getfreq(freqfn);
    p = loadfreq(freqfn);
    findsweeps(outfn, p, gridsize, 1, -1);
    free(p);
  }
  else if (argv[1][1]=='m') {
    int n, minx, maxx, x, j;
    if (argc!=5) usage();
    gridsize = atoi(argv[2]);
    if (gridsize <=0) {
      printf("gridsize should be > 0!\n"); usage();
    }
    sprintf(snpfn, argv[3]);
    sprintf(freqfn, "tempfreq.txt");
    sprintf(outfn, argv[4]);
    rep=0;
    while (readsnps_ms(snpfn)) {
      getfreq(freqfn);
      p = loadfreq(freqfn);
      findsweeps(outfn, p, gridsize, 0, ++rep);
      free(p);
      for (n=nmin; n<=nmax; n++) {
	if (invar==0) {minx=1; maxx=n;}
	else if (invar==1) {minx=1; maxx=n+1;}
	else {minx=0; maxx=n+1;}
	for (x=minx; x<maxx; x++) {
	  for (j=0; j<2; j++) 
	    free(prob[n][x][j]);
	  free(prob[n][x]);
	}
	free(prob[n]);
      }
      free(prob);
      prob=NULL;
    }
  }  
  else if (argv[1][1]=='M') {
    if (argc!=6) usage();
    gridsize = atoi(argv[2]);
    if (gridsize <= 0) {
      printf("gridsize should be > 0!\n"); usage();
    }
    sprintf(snpfn, argv[3]);
    sprintf(freqfn, argv[4]);
    sprintf(outfn, argv[5]);
    rep = 0;
    while (readsnps_ms(snpfn)) {
      p = loadfreq(freqfn);
      findsweeps(outfn, p, gridsize, 0, ++rep);
      free(p);
    }
  }
  else usage();
  return 0;
}
