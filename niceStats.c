/******* niceStats.c ********
for calculating sample stats from MS output 
********************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "miscCode/msGeneralStats.h"

#define LINEBUF 1000000

void usage();

int maxsites = 100000 ;

int main(argc,argv)
	int argc;
char *argv[];
{
	int nsam, i,  howmany  ,nwins;
	char **list, **cmatrix(), line[LINEBUF+1]  ;
	FILE *fopen(), *pfin ;
	double *posit ;
	int   segsites, count  , nadv, iss, h;
	double pi , th,  z, H, tajD, w, h1, h2, h12;
	char dum[50], astr[100] ;
        int *haplotype_counts;


/* read in first two lines of output  (parameters and seed) */
	pfin = stdin ;
	fgets( line, LINEBUF, pfin);
	sscanf(line," %s  %d %d", dum,  &nsam, &howmany);
	fgets( line, LINEBUF, pfin);

	if( argc > 1 ) { 
		nadv = atoi( argv[1] ) ; 
	}

	list = cmatrix(nsam,maxsites+1);
	posit = (double *)malloc( maxsites*sizeof( double ) ) ;

	count=0;
	//print header
	nwins=5;
	printf("pi\tss\tthetaH\ttajD\tfayWuH\tHapCount\tH1\tH12\tH2/H1\tOmega\tZnS");
	printf("\n");
	while( howmany-count++ ) {

		/* read in a sample */
		do {
			fgets( line, LINEBUF, pfin);
		}while ( line[0] != '/' );

		fscanf(pfin,"  segsites: %d", &segsites );
		if( segsites >= maxsites){
			maxsites = segsites + 10 ;
			posit = (double *)realloc( posit, maxsites*sizeof( double) ) ;
			biggerlist(nsam,maxsites, list) ;
		}
		if( segsites > 0) {
			fscanf(pfin," %s", astr);

			for( i=0; i<segsites ; i++) fscanf(pfin," %lf",posit+i) ;
			for( i=0; i<nsam;i++) fscanf(pfin," %s", list[i] );
		}
		/* analyse sample ( do stuff with segsites and list) */

	        haplotype_counts = (int *)malloc( nsam*sizeof( int ) ) ;

		iss = segSites(segsites,nsam,list);
		pi = nucdiv(nsam, segsites, list) ;
		th = thetah(nsam, segsites, list) ;
		h = nHaplotypes(segsites,nsam,list);
                getHaplotypeFreqSpec(segsites, nsam, list, haplotype_counts);
                h1 = petrovH1(haplotype_counts,nsam);
                h2 = petrovH2(haplotype_counts,nsam);
                h12 = petrovH12(haplotype_counts,nsam);
		H = th-pi;
		tajD = tajd(nsam,segsites,pi);
		z = ZnS( segsites,  nsam,  list);
                w = omegaMax(segsites, nsam,list);
		
		printf("%lf\t%d\t%lf\t%lf\t%lf\t%d\t%f\t%f\t%f\t%f\t%f", pi, iss, th ,tajD,H, h, h1, h12, h2/h1, w, z);
        free(haplotype_counts);
	printf("\n");
 }
	return(0);
}
