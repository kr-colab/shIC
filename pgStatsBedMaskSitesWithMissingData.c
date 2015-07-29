//pgStats for a bedFile rather than windows

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "miscCode/stringWrap.h"
#include "miscCode/sequenceMatrix.h"
#include "miscCode/pgSummaryStats.h"
#include "miscCode/bedFile.h"


void usage();

int main(int argc, char *argv[]){
	struct sequenceMatrix *aSeqMat, *ancSeqMat;
	int i, bedElNumber, h, ss;
	struct bedEl data[100000];
	double pi, theta_h, z, w, taj, H, h1, h2, h12;
        double maxFracMissingData;
	int numSites;
        int *haplotype_counts;
	
	if(argc < 5){
		usage();
		exit(1);
	}

	//open fastaFile and bedFile
	aSeqMat = sequenceMatrix_importFasta(argv[1]);
	ancSeqMat = sequenceMatrix_importFasta(argv[2]);
	bedElNumber = bedFileImport3(argv[3], data);
	maxFracMissingData = atof(argv[4]);

        haplotype_counts = (int *)malloc( aSeqMat->sampleSize*sizeof( int ) ) ;

	sequenceMatrix_NOutSitesWithTooMuchMissingData(aSeqMat,maxFracMissingData);

	//normal header
	printf("chrom\tchromStart\tchromEnd\tnumSites\tpi\tsegSites\tthetaH\ttajD\tfayWuH\tHapCount\tH1\tH12\tH2/H1\tOmega\tZnS\n");

	//loop through beds; the adjustments to end are to honor the zero indexed half open bed convention
	for(i=0;i<bedElNumber;i++){
                numSites = numColumnsNotNedOutFromTo(aSeqMat, data[i].chromStart, data[i].chromEnd);
		printf("%s\t%ld\t%ld\t%d\t",data[i].chrom,data[i].chromStart,data[i].chromEnd,numSites);
		pi = nucdivFromTo(aSeqMat, data[i].chromStart, data[i].chromEnd);
		ss = segSiteCountFromTo(aSeqMat, data[i].chromStart, data[i].chromEnd);
		taj = tajd(aSeqMat->sampleSize,ss,pi);
		theta_h = thetaHFromTo(aSeqMat,ancSeqMat,data[i].chromStart, data[i].chromEnd);
		H = theta_h - pi;
		//h = nHaplotypes(aSeqMat,data[i].chromStart, data[i].chromEnd);
                h = getHaplotypeFreqSpec(aSeqMat, data[i].chromStart, data[i].chromEnd, haplotype_counts);
                h1 = petrovH1(haplotype_counts, aSeqMat->sampleSize);
                h2 = petrovH2(haplotype_counts, aSeqMat->sampleSize);
                h12 = petrovH12(haplotype_counts, aSeqMat->sampleSize);
		z = ZnSFromTo(aSeqMat,data[i].chromStart, data[i].chromEnd);
		w = omegaMaxFromTo(aSeqMat,data[i].chromStart,data[i].chromEnd);
		printf("%f\t%d\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\n",pi, ss, theta_h, taj, H, h, h1, h12, h2/h1, w, z);
	}
	sequenceMatrix_free(aSeqMat);
	return(0);
}	

void usage(){
	printf("pgStatsBed ingroupFastaFile ancestorFastaFile bedFile maxFractionMissingData\n");
}
