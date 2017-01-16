/***** ms.c     ************************************************
*
*       Generates samples of gametes ( theta given or fixed number 
*						of segregating sites.)
*	Usage is shown by typing ms without arguments.   
        usage: ms nsam howmany  -t  theta  [options]
		or
	       ms nsam howmany -s segsites  [options] 

	   nsam is the number of gametes per sample.
	   howmany is the number of samples to produce.
	   With -t the numbers of segregating sites will randomly vary 
		from one sample to the next.
	   with -s segsites,  the number of segregating sites will be
		segsites in each sample.

           Other options: See msdoc.pdf or after downloading and compiling, type ms<CR>.


*	  Arguments of the options are explained here:

           npop:  Number of subpopulations which make up the total population
           ni:  the sample size from the i th subpopulation (all must be 
		specified.) The output will have the gametes in order such that
		the first n1 gametes are from the first island, the next n2 are
		from the second island, etc.
           nsites: number of sites between which recombination can occur.
           theta: 4No times the neutral mutation rate 
           rho: recombination rate between ends of segment times 4No
	   f: ratio of conversion rate to recombination rate. (Wiuf and Hein model.)
	   track_len:  mean length of conversion track in units of sites.  The 
		       total number of sites is nsites, specified with the -r option.
           mig_rate: migration rate: the fraction of each subpop made up of
                 migrants times 4No. 
           howmany: howmany samples to generate.

	Note:  In the above definition, No is the total diploid population if
		npop is one, otherwise, No is the diploid population size of each
		subpopulation. 
	A seed file called "seedms" will be created  if it doesn't exist. The
		seed(s) in this file will be modified by the program. 
		So subsequent runs
		will produce new output.  The initial contents of seedms will be
		printed on the second line of the output.
        Output consists of one line with the command line arguments and one
	 	line with the seed(s).
		The samples appear sequentially following that line.
		Each sample begins with "//", then the number of segregating sites, the positions
		of the segregating sites (on a scale of 0.0 - 1.0). On the following
		lines are the sampled gametes, with mutants alleles represented as
		ones and ancestral alleles as zeros.
	To compile:  cc -o ms  ms.c  streec.c  rand1.c -lm
		or:  cc -o ms ms.c streec.c rand2.c -lm
	 (Of course, gcc would be used instead of cc on some machines.  And -O3 or 
		some other optimization switches might be usefully employed with some 
		compilers.) ( rand1.c uses drand48(), whereas rand2.c uses rand() ).

*
*   Modifications made to combine ms and mss on 25 Feb 2001
*	Modifications to command line options to use switches  25 Feb 2001
*	Modifications to add // before each sample  25 Feb 2001
	Modifications to add gene conversion 5 Mar 2001
	Added demographic options -d  13 Mar 2001
	Changed ran1() to use rand(). Changed seed i/o to accomodate this change. 20 April.
	Changed cleftr() to check for zero rand() .13 June 2001
	Move seed stuff to subroutine seedit()  11 July 2001
	Modified streec.c to handle zero length demographic intervals 9 Aug 2001
	Corrected problem with negative growth rates (Thanks to D. Posada and C. Wiuf) 13 May 2002
	Changed sample_stats.c to output thetah - pi rather than pi - thetah.  March 8 2003.
	Changed many command line options, allowing arbitrary migration matrix, and subpopulation
	   sizes.  Also allows parameters to come from a file. Option to output trees.  Option to
	   split and join subpopulations.   March 8, 2003. (Old versions saved in msold.tar ).
	!!! Fixed bug in -en option.  Earlier versions may have produced garbage when -en ... used. 9 Dec 2003
	Fixed bug which resulted in incorrect results for the case where
             rho = 0.0 and gene conversion rate > 0.0. This case was not handled
	    correctly in early versions of the program. 5 Apr 2004.  (Thanks to
	    Vincent Plagnol for pointing out this problem.) 
	Fixed bug in prtree().  Earlier versions may have produced garbage when the -T option was used.
		 1 Jul 2004.
	Fixed bug in -e. options that caused problems with -f option  13 Aug 2004.
	Fixed bug in -es option, which was a problem when used with -eG. (Thanks again to V. Plagnol.) 6 Nov. 2004
	Added -F option:  -F minfreq  produces output with sites with minor allele freq < minfreq filtered out.  11 Nov. 2004.
	Fixed bug in streec.c (isseg() ).  Bug caused segmentation fault, crash on some machines. (Thanks
	    to Melissa Jane Todd-Hubisz for finding and reporting this bug.)
***************************************************************************/


#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "ms.h"

#define SITESINC 10 
#define DEBUG 0
#define DEBUG1 1

unsigned maxsites = SITESINC ;

struct node{
	int abv;
	int ndes;
	float time;
	};

struct segl {
	int beg;
	struct node *ptree;
	int next;
	};

double *posit ;
double segfac ;
	
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
	
	bool *out1;
    double *out2;
    int nsites, segs, nsam, i, j, howmany, npop, *config, pop, count, narg; 
    double  r,theta, mig_rate, estimator, gst, alpha, within, between;
    char **list, **cmatrix();
    double ttot, f, track_len, probss;
    int  segsites, segsitesin, nintn;
	double *nrec, *npast, *tpast ;
	
	void seedit( const char * ) ;
	struct params param, getdefpars( int nsam, double theta, int howmany, double r, int nsites ) ;
	int gensam( struct params *param, char **list, double *probss ) ;
	struct devent *pevent ;    
    
  if(nrhs<2) {
	  mexErrMsgTxt("At least two inputs required, e.g.,MS_MEX(10,4.0)\n[gametes,position]=MS_MEX(nsam,theta,segs,r,nsites)");
  } else if(nlhs>2) {
     mexErrMsgTxt("Too many output arguments");
  }
  
  /* NSAM */
  if (mxIsEmpty(prhs[0])|| (!mxIsNumeric(prhs[0]))){
    mexErrMsgTxt("NSAM must be an integer.");
  }
  nsam = (int)mxGetScalar(prhs[0]);
  if(nsam<2) {
    mexErrMsgTxt("nsam must >=2.");
  }

  /* THETA */
  if (mxIsEmpty(prhs[1])|| (!mxIsNumeric(prhs[1]))){
    mexErrMsgTxt("THETA must be a double.");
  }
  theta = mxGetScalar(prhs[1]);

  /* SEGS */
  if (nrhs>2){
	  if (mxIsEmpty(prhs[2])|| (!mxIsNumeric(prhs[2]))){
		mexErrMsgTxt("SEGS must be an integer.");
	  }
	segs = (int)mxGetScalar(prhs[2]);
  }else{
	segs=0;
  }

if (segs==0 && theta==0){
		mexErrMsgTxt("THETA and SEGS cannot be zero at same time.");
}

  /* R */
  if (nrhs>3){
	  if (mxIsEmpty(prhs[3])|| (!mxIsNumeric(prhs[3]))){
		mexErrMsgTxt("R must be a double.");
	  }
	r = mxGetScalar(prhs[3]);
  }else{
	r=0.0;
  }

  /* NSITES */ 
  if (nrhs>4){
	  if (mxIsEmpty(prhs[4])|| (!mxIsNumeric(prhs[4]))){
		mexErrMsgTxt("NSITES must be an integer.");
	  }
	nsites = (int)mxGetScalar(prhs[4]);
  }else{
	nsites=2;
  }

  param = getdefpars(nsam, theta, segs, r, nsites) ;

if( (param.mp.theta == 0.0) && ( param.mp.segsitesin == 0 ) && ( param.mp.treeflag == 0 ) ) {
	mexErrMsgTxt(" either -s or -t or -T option must be used.");
}

/* seedit("s"); */

/* howmany=1; pf = stdout ; */

    if( param.mp.segsitesin ==  0 ) {
	   list = cmatrix(param.cp.nsam,maxsites+1);
	   posit = (double *)malloc( (unsigned)( maxsites*sizeof( double)) ) ;
	   }
	else {
	    list = cmatrix(param.cp.nsam, param.mp.segsitesin+1 ) ;
	    posit = (double *)malloc( (unsigned)( param.mp.segsitesin*sizeof( double)) ) ;
	    if( param.mp.theta > 0.0 ){
		segfac = 1.0 ;
		for(  i= param.mp.segsitesin; i > 1; i--) segfac *= i ;
		}
	   }


  segsites = gensam( &param, list, &probss ) ; 
  plhs[0] = mxCreateLogicalMatrix(nsam,segsites);
  plhs[1] = mxCreateDoubleMatrix(1,segsites,mxREAL);
   
  out2  = mxGetPr(plhs[1]);
  for( i=0; i<segsites; i++){
			if(DEBUG){
				mexPrintf("%6.4lf ",posit[i] ); mexPrintf("\n");
			}
            out2[i]= posit[i];
			}
  
  out1 = (bool *)mxGetLogicals(plhs[0]);
  
		   for( i=0; i<nsam; i++){
           for( j=0; j<segsites; j++){
			 count=i+nsam*j;
			 if ((list[i][j]-48)==1){
				if(DEBUG) mexPrintf("(%d,%d)=%d==%d\n",i+1,j+1,count+1,list[i][j]-48);
				out1[count]=true;
			 }
		   }}
		   
		   /*  seedit( "end" ); */
           return;
}


int 
gensam( struct params *pparam, char **list, double *pprobss ) 
{
	int nsegs, h, i, k, j, seg, ns, start, end, len, segsit ;
	struct segl *seglst, *segtre_mig(struct c_params *p, int *nsegs ) ; /* used to be: [MAXSEG];  */
	double nsinv,  tseg, tt, ttime(struct node *, int nsam), ttimemf(struct node *, int nsam, int mfreq) ;
	double *pk;
	int *ss;
	int segsitesin,nsites;
	double theta, es ;
	int nsam, mfreq ;
	void prtree( struct node *ptree, int nsam);
        int make_gametes(int nsam, int mfreq,  struct node *ptree, double tt, int newsites, int ns, char **list );
 	void ndes_setup( struct node *, int nsam );


	nsites = pparam->cp.nsites ;
	nsinv = 1./nsites;
	seglst = segtre_mig(&(pparam->cp),  &nsegs ) ;
	
	nsam = pparam->cp.nsam;
	segsitesin = pparam->mp.segsitesin ;
	theta = pparam->mp.theta ;
	mfreq = pparam->mp.mfreq ;

	if( pparam->mp.treeflag ) {
	  	ns = 0 ;
	    for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
	      if( (pparam->cp.r > 0.0 ) || (pparam->cp.f > 0.0) ){
		end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
		start = seglst[seg].beg ;
		len = end - start + 1 ;
	/*	fprintf(stdout,"[%d]", len); */
	      }
	      prtree( seglst[seg].ptree, nsam ) ;
	      if( (segsitesin == 0) && ( theta == 0.0 ) ) 
	  	    free(seglst[seg].ptree) ;
	    }
	}

    if( (segsitesin == 0) && ( theta > 0.0)   ) {
	ns = 0 ;
	for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) { 
		if( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam );
		end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
		start = seglst[seg].beg ;
		len = end - start + 1 ;
		tseg = len*(theta/nsites) ;
		if( mfreq == 1) tt = ttime(seglst[seg].ptree, nsam);
                else tt = ttimemf(seglst[seg].ptree, nsam, mfreq );
		segsit = poisso( tseg*tt );
		if( (segsit + ns) >= maxsites ) {
			maxsites = segsit + ns + SITESINC ;
			posit = (double *)realloc(posit, maxsites*sizeof(double) ) ;
			  biggerlist(nsam, list) ; 
			}
		make_gametes(nsam,mfreq,seglst[seg].ptree,tt, segsit, ns, list );
		free(seglst[seg].ptree) ;
		locate(segsit,start*nsinv, len*nsinv,posit+ns);   
		ns += segsit;
		}
    }
   else if( segsitesin > 0 ) {

        pk = (double *)malloc((unsigned)(nsegs*sizeof(double)));
        ss = (int *)malloc((unsigned)(nsegs*sizeof(int)));
        if( (pk==NULL) || (ss==NULL) ) perror("malloc error. gensam.2");


	tt = 0.0 ;
	for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) { 
		if( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam );
		end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
		start = seglst[seg].beg ;
		len = end - start + 1 ;
		tseg = len/(double)nsites ;
               if( mfreq == 1 ) pk[k] = ttime(seglst[seg].ptree,nsam)*tseg ;
               else pk[k] = ttimemf(seglst[seg].ptree,nsam, mfreq)*tseg ;
                 tt += pk[k] ;
		}
	if( theta > 0.0 ) { 
	  es = theta * tt ;
	  *pprobss = exp( -es )*pow( es, (double) segsitesin) / segfac ;
	  }
       if( tt > 0.0 ) {
          for (k=0;k<nsegs;k++) pk[k] /= tt ;
          mnmial(segsitesin,nsegs,pk,ss);
          }
        else
            for( k=0; k<nsegs; k++) ss[k] = 0 ;
	ns = 0 ;
	for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) { 
		end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
		start = seglst[seg].beg ;
		len = end - start + 1 ;
                tseg = len/(double)nsites;
                make_gametes(nsam,mfreq,seglst[seg].ptree,tt*pk[k]/tseg, ss[k], ns, list);

		free(seglst[seg].ptree) ;
		locate(ss[k],start*nsinv, len*nsinv,posit+ns);   
		ns += ss[k] ;
		}
	free(pk);
	free(ss);

    }
	for(i=0;i<nsam;i++) list[i][ns] = '\0' ;
	return( ns ) ;
}

	void 
ndes_setup(struct node *ptree, int nsam )
{
	int i ;

	for( i=0; i<nsam; i++) (ptree+i)->ndes = 1 ;
	for( i = nsam; i< 2*nsam -1; i++) (ptree+i)->ndes = 0 ;
	for( i= 0; i< 2*nsam -2 ; i++)  (ptree+((ptree+i)->abv))->ndes += (ptree+i)->ndes ;

}

	int
biggerlist(nsam,  list )
	int nsam ;
	char ** list ;
{
	int i;

/*  fprintf(stderr,"maxsites: %d\n",maxsites);  */	
	for( i=0; i<nsam; i++){
	   list[i] = (char *)realloc( list[i],maxsites*sizeof(char) ) ;
	   if( list[i] == NULL ) perror( "realloc error. bigger");
	   }
}
	   


/* allocates space for gametes (character strings) */
	char **
cmatrix(nsam,len)
	int nsam, len;
{
	int i;
	char **m;

	if( ! ( m = (char **) malloc( (unsigned) nsam*sizeof( char* ) ) ) )
	   perror("alloc error in cmatrix") ;
	for( i=0; i<nsam; i++) {
		if( ! ( m[i] = (char *) malloc( (unsigned) len*sizeof( char ) )))
			perror("alloc error in cmatric. 2");
		}
	return( m );
}



	int
locate(n,beg,len,ptr)
	int n;
	double beg, len, *ptr;
{
	int i;

	ordran(n,ptr);
	for(i=0; i<n; i++)
		ptr[i] = beg + ptr[i]*len ;

}



  struct params
getdefpars(int nsam, double theta, int segs, double r, int nsites)
{
	int arg, i, j, sum , pop , argstart, npop , npop2, pop2 ;
	double migr, mij, psize, palpha ;
	struct params p;
	struct devent *ptemp , *pt ;
	arg = 3 ;
	p.cp.r = p.mp.theta =  p.cp.f = 0.0 ;

	p.mp.theta = theta;
	p.cp.nsam = nsam;
    


	p.cp.track_len = 0. ;
	p.cp.npop = npop = 1 ;
	p.cp.mig_mat = (double **)malloc( (unsigned) sizeof( double *) );
	p.cp.mig_mat[0] = (double *)malloc( (unsigned)sizeof(double *));
	p.cp.mig_mat[0][0] =  0.0 ;
	p.mp.segsitesin = 0 ;

    p.mp.segsitesin = (int)segs;

	p.mp.treeflag = 0 ;
        p.mp.mfreq = 1 ;
	p.cp.config = (int *) malloc( (unsigned)(( p.cp.npop +1 ) *sizeof( int)) );
	(p.cp.config)[0] = p.cp.nsam ;
	p.cp.size= (double *) malloc( (unsigned)( p.cp.npop *sizeof( double )) );
	(p.cp.size)[0] = 1.0  ;
	p.cp.alphag = (double *) malloc( (unsigned)(( p.cp.npop ) *sizeof( double )) );
	(p.cp.alphag)[0] = 0.0  ;
	p.cp.r=r;
	if (r==0){
	  p.cp.nsites = 2 ;
	}else{
	  p.cp.nsites = (int)nsites;
	}

	p.cp.deventlist = NULL ;

		if( (p.mp.theta == 0.0) && ( p.mp.segsitesin == 0 ) && ( p.mp.treeflag == 0 ) ) {
			fprintf(stderr," either -s or -t or -T option must be used. \n");
			usage();
			exit(1);
			}
		sum = 0 ;
		for( i=0; i< p.cp.npop; i++) sum += (p.cp.config)[i] ;
		if( sum != p.cp.nsam ) {
			fprintf(stderr," sum sample sizes != nsam\n");
			usage();
			exit(1);
			}


	return p ;
}


	void
argcheck( int arg, int argc, char *argv[] )
{
	if( (arg >= argc ) || ( argv[arg][0] == '-') ) {
	   fprintf(stderr,"not enough arguments after %s\n", argv[arg-1] ) ;
	   fprintf(stderr,"For usage type: ms<return>\n");
	   exit(0);
	  }
}
	
	int
usage()
{
fprintf(stderr,"usage: ms nsam howmany \n");
fprintf(stderr,"  Options: \n"); 
fprintf(stderr,"\t -t theta   (this option and/or the next must be used. Theta = 4*N0*u )\n");
fprintf(stderr,"\t -s segsites   ( fixed number of segregating sites)\n");
fprintf(stderr,"\t -T          (Output gene tree.)\n");
fprintf(stderr,"\t -F minfreq     Output only sites with freq of minor allele >= minfreq.\n");
fprintf(stderr,"\t -r rho nsites     (rho here is 4Nc)\n");
fprintf(stderr,"\t\t -c f track_len   (f = ratio of conversion rate to rec rate. tracklen is mean length.) \n");
fprintf(stderr,"\t\t\t if rho = 0.,  f = 4*N0*g, with g the gene conversion rate.\n"); 
fprintf(stderr,"\t -G alpha  ( N(t) = N0*exp(-alpha*t) .  alpha = -log(Np/Nr)/t\n");      
fprintf(stderr,"\t -I npop n1 n2 ... [mig_rate] (all elements of mig matrix set to mig_rate/(npop-1) \n");    
fprintf(stderr,"\t\t -m i j m_ij    (i,j-th element of mig matrix set to m_ij.)\n"); 
fprintf(stderr,"\t\t -ma m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n"); 
fprintf(stderr,"\t\t -n i size_i   (popi has size set to size_i*N0 \n");
fprintf(stderr,"\t\t -g i alpha_i  (If used must appear after -M option.)\n"); 
fprintf(stderr,"\t   The following options modify parameters at the time 't' specified as the first argument:\n");
fprintf(stderr,"\t -eG t alpha  (Modify growth rate of all pop's.)\n");     
fprintf(stderr,"\t -eg t i alpha_i  (Modify growth rate of pop i.) \n");    
fprintf(stderr,"\t -eM t mig_rate   (Modify the mig matrix so all elements are mig_rate/(npop-1)\n"); 
fprintf(stderr,"\t -em t i j m_ij    (i,j-th element of mig matrix set to m_ij at time t )\n"); 
fprintf(stderr,"\t -ema t npop  m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n");  
fprintf(stderr,"\t -eN t size  (Modify pop sizes. New sizes = size*N0 ) \n");    
fprintf(stderr,"\t -en t i size_i  (Modify pop size of pop i.  New size of popi = size_i*N0 .)\n");
fprintf(stderr,"\t -es t i proportion  (Split: pop i -> pop-i + pop-npop, npop increases by 1.\n");    
fprintf(stderr,"\t\t proportion is probability that each lineage stays in pop-i. (p, 1-p are admixt. proport.\n");
fprintf(stderr,"\t\t Size of pop npop is set to N0 and alpha = 0.0 , size and alpha of pop i are unchanged.\n");
fprintf(stderr,"\t -ej t i j   ( Join lineages in pop i and pop j into pop j\n");
fprintf(stderr,"\t\t  size, alpha and M are unchanged.\n");  
fprintf(stderr,"\t  -f filename     ( Read command line arguments from file filename.)\n");  
fprintf(stderr," See msdoc.pdf for explanation of these parameters.\n");

exit(1);
}
	void
addtoelist( struct devent *pt, struct devent *elist ) 
{
	struct devent *plast, *pevent, *ptemp  ;

	pevent = elist ;
	while(  (pevent != NULL ) && ( pevent->time <= pt->time ) )  {
		plast = pevent ;
		pevent = pevent->nextde ;
		}
	ptemp = plast->nextde ;
	plast->nextde = pt ;
	pt->nextde = ptemp ;
}
	
/************ make_gametes.c  *******************************************
*
*
*****************************************************************************/

#define STATE1 '1'
#define STATE2 '0'

	int
make_gametes(nsam,mfreq, ptree,tt,newsites, ns, list )
	int nsam,mfreq, newsites, ns;
	struct node *ptree;
	char **list;
	double tt ;
{
	int  tip, j,  node ;
        int pickb(int nsam, struct node *ptree, double tt), 
            pickbmf(int nsam, int mfreq, struct node *ptree, double tt) ;

	for(  j=ns; j< ns+newsites ;  j++ ) {
		if( mfreq == 1 ) node = pickb(  nsam, ptree, tt);
		else node = pickbmf(  nsam, mfreq, ptree, tt);
		for( tip=0; tip < nsam ; tip++) {
		   if( tdesn(ptree, tip, node) ) list[tip][j] = STATE1 ;
		   else list[tip][j] = STATE2 ;
		   }
		}
}


/***  ttime.c : Returns the total time in the tree, *ptree, with nsam tips. **/

	double
ttime( ptree, nsam)
	struct node *ptree;
	int nsam;
{
	double t;
	int i;

	t = (ptree + 2*nsam-2) -> time ;
	for( i=nsam; i< 2*nsam-1 ; i++)
		t += (ptree + i)-> time ;
	return(t);
}


	double
ttimemf( ptree, nsam, mfreq)
	struct node *ptree;
	int nsam, mfreq;
{
	double t;
	int i;

	t = 0. ;
	for( i=0;  i< 2*nsam-2  ; i++)
	  if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) )
		t += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
	return(t);
}


	void
prtree( ptree, nsam)
	struct node *ptree;
	int nsam;
{
	double t;
	int i, *descl, *descr ;
	void parens( struct node *ptree, int *descl, int *descr, int noden );

	descl = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );
	descr = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );
	for( i=0; i<2*nsam-1; i++) descl[i] = descr[i] = -1 ;
	for( i = 0; i< 2*nsam-2; i++){
	  if( descl[ (ptree+i)->abv ] == -1 ) descl[(ptree+i)->abv] = i ;
	  else descr[ (ptree+i)->abv] = i ;
	 }
	parens( ptree, descl, descr, 2*nsam-2);
	free( descl ) ;
	free( descr ) ;
}

	void
parens( struct node *ptree, int *descl, int *descr,  int noden)
{
	double time ;

   if( descl[noden] == -1 ) {
	printf("%d:%5.3lf", noden+1, (ptree+ ((ptree+noden)->abv))->time );
	}
   else{
	printf("(");
	parens( ptree, descl,descr, descl[noden] ) ;
	printf(",");
	parens(ptree, descl, descr, descr[noden] ) ;
	if( (ptree+noden)->abv == 0 ) printf(");\n"); 
	else {
	  time = (ptree + (ptree+noden)->abv )->time - (ptree+noden)->time ;
	  printf("):%5.3lf", time );
	  }
        }
}

/***  pickb : returns a random branch from the tree. The probability of picking
              a particular branch is proportional to its duration. tt is total
	      time in tree.   ****/

	int
pickb(nsam, ptree, tt)
	int nsam;
	struct node *ptree;
	double tt;
{
	double x, y, ran1();
	int i;

	x = ran1()*tt;    
	for( i=0, y=0; i < 2*nsam-2 ; i++) {
		y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
		if( y >= x ) return( i ) ;
		}
	return( i );
}

	int
pickbmf(nsam, mfreq, ptree, tt )
	int nsam, mfreq;
	struct node *ptree;
	double tt;
{
	double x, y, ran1();
	int i;

	x = ran1()*tt;
	for( i=0, y=0; i < 2*nsam-2 ; i++) {
	  if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) )
		y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
	  if( y >= x ) return( i ) ;
	}
	return( i );
}

/****  tdesn : returns 1 if tip is a descendant of node in *ptree, otherwise 0. **/

	int
tdesn(ptree, tip, node )
	struct node *ptree;
	int tip, node;
{
	int k;

	for( k= tip ; k < node ; k = (ptree+k)->abv ) ;
	if( k==node ) return(1);
	else return(0);
}


/* pick2()  */

	int
pick2(n,i,j)
	int n, *i, *j;
{
	double ran1();

	*i = n * ran1() ;
	while( ( *j = n * ran1() ) == *i )
		;
	return(0) ;
}

/**** ordran.c  ***/

	int
ordran(n,pbuf)
	int n;
	double pbuf[];
{
	ranvec(n,pbuf);
	order(n,pbuf);
	/* return; */
}


	int
mnmial(n,nclass,p,rv)
	int n, nclass, rv[];
	double p[];
{
	double ran1();
	double x, s;
	int i, j;

	for(i=0; i<nclass; i++) rv[i]=0;
	for(i=0; i<n ; i++) {
	   x = ran1();
	   j=0;
	   s = p[0];
	   while( (x > s) && ( j<(nclass-1) ) )  s += p[++j];
	   rv[j]++;
	   }
	return(j);
}

        int
order(n,pbuf)
        int n;
        double pbuf[];
{
        int gap, i, j;
        double temp;

        for( gap= n/2; gap>0; gap /= 2)
           for( i=gap; i<n; i++)
                for( j=i-gap; j>=0 && pbuf[j]>pbuf[j+gap]; j -=gap) {
                   temp = pbuf[j];
                   pbuf[j] = pbuf[j+gap];
                   pbuf[j+gap] = temp;
                   }
}


	int
ranvec(n,pbuf)
	int n;
	double pbuf[];
{
	int i;
	double ran1();

	for(i=0; i<n; i++)
		pbuf[i] = ran1();

	/* return; */
}



	int
poisso(u)
	double u;
{
	double  cump, ru, ran1(), p, gasdev() ;
	int i=1;

	if( u > 30. ) return( (int)(0.5 + gasdev(u,u)) );
	ru = ran1();
	p = exp(-u);
	if( ru < p) return(0);
	cump = p;
	
	while( ru > ( cump += (p *= u/i ) ) )
		i++;
	return(i);
}


/* a slight modification of crecipes version */

double gasdev(m,v)
	double m, v;
{
	static int iset=0;
	static float gset;
	float fac,r,v1,v2;
	double ran1();

	if  (iset == 0) {
		do {
			v1=2.0*ran1()-1.0;
			v2=2.0*ran1()-1.0;
			r=v1*v1+v2*v2;
		} while (r >= 1.0);
		fac=sqrt(-2.0*log(r)/r);
		gset= v1*fac;
		iset=1;
		return( m + sqrt(v)*v2*fac);
	} else {
		iset=0;
		return( m + sqrt(v)*gset ) ;
	}
}

