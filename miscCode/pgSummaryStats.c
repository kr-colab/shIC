/*****************************
/ pgSummaryStats.c 
/ 
/ A. D. Kern
/
/ popGen statistics for
/ sequenceMatrix objects
******************************/


#include "pgSummaryStats.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "vector.h"

#define MAXMUTS 1000

int nHaplotypes(struct sequenceMatrix* aSeqMat, int beg, int end)
{
	int i;
	int j;
	int k;
	int haplotype_found;
	int allsame;
	
	int n_haplotypes = 0;
	char **haplotypes; //haplotypes[aSeqMat->sampleSize][end-beg];
	
	 if( ! ( haplotypes = (char **) malloc( (unsigned)( aSeqMat->sampleSize*sizeof( char* )) ) ) )
		perror("alloc error in haplotypes") ;
	for( i=0; i<aSeqMat->sampleSize; i++) {
		if( ! ( haplotypes[i] = (char *) malloc( (unsigned) ((end-beg)*sizeof( char )) )))
			perror("alloc error in haplotypes 2");
	}

	for(i=0; i<aSeqMat->sampleSize; i++)
	{
		haplotype_found = 0;
		for(j=0; j<n_haplotypes; j++)
		{
			allsame = 1;
			for(k=0; k<end-beg; k++)
			{
				if(haplotypes[j][k] != aSeqMat->matrix[i]->cString[beg+k])
				{
					if(haplotypes[j][k] == 'N' || haplotypes[j][k] == '-')
						haplotypes[j][k] = aSeqMat->matrix[i]->cString[beg+k];
												
					else if(aSeqMat->matrix[i]->cString[beg+k] != 'N' && aSeqMat->matrix[i]->cString[beg+k] != '-')
					{
						allsame = 0;
						break;	
					}
				}
			}
			
			if(allsame)
			{
				haplotype_found = 1;
				break;	
			}
		}
		
		if(!haplotype_found)
		{
			n_haplotypes++;
			for(j=0; j<end-beg; j++)
					haplotypes[n_haplotypes-1][j] = aSeqMat->matrix[i]->cString[beg+j];
		}
	}
	for(i=0;i<aSeqMat->sampleSize;i++){
		free(haplotypes[i]);
	}	
	free(haplotypes);
	return n_haplotypes;
}

//haplotype_counts[i] is the number of haplotypes found in exactly i+1 individuals
int getHaplotypeFreqSpec(struct sequenceMatrix* aSeqMat, int beg, int end, int *haplotype_counts)
{
        int i;
        int j;
        int k;
        int haplotype_found;
        int allsame;
        int freq;

        int n_haplotypes = 0;
        int haplotype_occurrences[aSeqMat->sampleSize];
        char **haplotypes; //haplotypes[aSeqMat->sampleSize][end-beg];
        
        if( ! ( haplotypes = (char **) malloc( (unsigned)( aSeqMat->sampleSize*sizeof( char* )) ) ) )
                perror("alloc error in haplotypes") ;
        for( i=0; i<aSeqMat->sampleSize; i++) {
                if( ! ( haplotypes[i] = (char *) malloc( (unsigned) (((end-beg)+1)*sizeof( char )) )))
                        perror("alloc error in haplotypes 2");
        }

        for(i=0; i<aSeqMat->sampleSize; i++)
        {
            haplotype_counts[i] = 0;
            haplotype_occurrences[i] = 0;
        }

        for(i=0; i<aSeqMat->sampleSize; i++)
        {
                haplotype_found = 0;
                for(j=0; j<n_haplotypes; j++)
                {
                        allsame = 1;
                        for(k=0; k<end-beg; k++)
                        {
                                if(haplotypes[j][k] != aSeqMat->matrix[i]->cString[beg+k])
                                {
                                        if(haplotypes[j][k] == 'N' || haplotypes[j][k] == '-')
                                                haplotypes[j][k] = aSeqMat->matrix[i]->cString[beg+k];
                                                                                                
                                        else if(aSeqMat->matrix[i]->cString[beg+k] != 'N' && aSeqMat->matrix[i]->cString[beg+k] != '-')
                                        {
                                                allsame = 0;
                                                break;  
                                        }
                                }
                        }

                        if(allsame)
                        {
                                haplotype_found = 1;
                                haplotype_occurrences[j]+=1;
                                break;
                        }
                }
                if(!haplotype_found)
                {
                        n_haplotypes++;
                        for(j=0; j<end-beg; j++)
                                        haplotypes[n_haplotypes-1][j] = aSeqMat->matrix[i]->cString[beg+j];
                        haplotypes[n_haplotypes-1][end-beg]='\0';
                        haplotype_occurrences[n_haplotypes-1]=1;
                }
        }

        for (i=0; i<n_haplotypes; i++)
        {
                freq = haplotype_occurrences[i];
                if (freq > 0 && freq <= aSeqMat->sampleSize)
                {
                        haplotype_counts[freq-1] += 1;
                }
        }

        for(i=0;i<aSeqMat->sampleSize;i++){
                free(haplotypes[i]);
        }       
        free(haplotypes);
        return n_haplotypes;
}

double petrovH1(int *haplotype_counts, int nsam)
{
    int hapFreq;
    double pi;
    double h1 = 0.0;

    for (hapFreq=nsam; hapFreq>0; hapFreq--)
    {
        pi = hapFreq/ (double)nsam;
        h1 += haplotype_counts[hapFreq-1]*pi*pi;
    }
    return h1;
}

double petrovH2(int *haplotype_counts, int nsam)
{
    int hapFreq;
    double pi;
    double h2 = 0.0;
    int first = 1;

    for (hapFreq=nsam; hapFreq>0; hapFreq--)
    {
        pi = hapFreq/ (double)nsam;
        if (haplotype_counts[hapFreq-1] > 0)
        {
            if (first)
            {
                first = 0;
                h2 += (haplotype_counts[hapFreq-1]-1)*pi*pi;
            }
            else
            {
                h2 += haplotype_counts[hapFreq-1]*pi*pi;
            }
        }
    }
    return h2;
}

double petrovH12(int *haplotype_counts, int nsam)
{
    int hapFreq, i;
    double pi;
    double part1 = 0.0;
    double part2 = 0.0;
    int totalAdded = 0;

    for (hapFreq=nsam; hapFreq>0; hapFreq--)
    {
        pi = hapFreq/ (double)nsam;
        for (i = 0;i < haplotype_counts[hapFreq-1];i++)
        {
            if (totalAdded < 2)
            {
                part1 += pi;
            }
            else
            {
                part2 += pi*pi;

            }
            totalAdded++;
        }
    }

    part1 = part1*part1;
    return part1+part2;
}

double nucdiv(struct sequenceMatrix *aSeqMat){
	double pi = 0.0;
	pi = nucdivFromTo(aSeqMat,0,aSeqMat->length);
	return(pi);	
}

// fixed differences (set is ingroup data)
int fixedDiffsPre(stringWrap** set, int size, int samplesize, stringWrap* outgdata)
{
	int i;
	int fixeddiffs = 0;

	for(i=0; i<size; i++)
	{
		if(set[i]->length == 1)
		{
			/* if outgdata != set[0] it's a fixed difference */
			if(outgdata->cString[0] != '-' && outgdata->cString[0] != 'N' && outgdata->cString[0] != set[i]->cString[0])
				fixeddiffs++;
		}
	}

	return fixeddiffs;
}

// proportion of missing data
double missingDataPre(stringWrap** array, int size, int samplesize)
{
	int i;
	int total;
	int present;

	total = samplesize * size;
	present = 0;

	for(i=0; i<size; i++)
		present += array[i]->length;

	return (double)(total - present) / (double)(total);
}

double missingDataFromTo(struct sequenceMatrix* aSeqMat, int start, int end)
{
	int i;
	stringWrap* array = stringWrapNew(aSeqMat->sampleSize);

	int total = (end-start) * aSeqMat->sampleSize;
	int present = 0;

	for(i=start; i<end; i++)
	{
		sequenceMatrix_siteArrayClean(aSeqMat, array, i);
		present += array->length;
	}

	stringWrapFree(array);
	return (double)(total - present) / (double)total;
}

/* calculates pi (nucdiv) using pre-calculated set and array */
double nucdivPre(stringWrap** set, stringWrap** array, int size)
{
	int i, j;
	double pi, p1, n, nnm1, sum;

	pi  = 0.0;
	sum = 0.0;

	for(i=0; i<size; i++)
	{
		if(set[i]->length > 1)
		{
			n = (double)(array[i]->length);
			nnm1 = n / (n-1.0);
			p1 = 0.0;
			pi = 0.0;

			for(j=0; j<set[i]->length; j++)
			{
				p1 = (double)(stringWrapCountChar(array[i], set[i]->cString[j])) / n;
				pi += p1*p1;
			}
			
			sum += (1.0-pi) * nnm1;
		}
	}

	return sum;
}

//redundant with piWindow below but easier to use
double nucdivFromTo(struct sequenceMatrix *aSeqMat, int start, int end){
	int i, j;
  	double pi, p1, n, nnm1,sum ;
	stringWrap *set, *array;

  	pi = 0.0 ;
	sum =0.0;
	set = stringWrapNew(aSeqMat->sampleSize);
	array = stringWrapNew(aSeqMat->sampleSize);
  	for( i = start; i < end; i++){
		sequenceMatrix_siteSetClean(aSeqMat,set,i);
		if(set->length > 1){
			sequenceMatrix_siteArrayClean(aSeqMat,array,i);
			n = (double) array->length;
			nnm1 = n / (n-1.0);
			p1 = 0.;
			pi = 0.;
			for(j=0; j < set->length; j++){
				p1 = (double) stringWrapCountChar(array,set->cString[j]) / n;
				pi += (p1*p1);
			}
			sum += (1.0-pi) * nnm1;
  		}
	}
	stringWrapFree(set);
	stringWrapFree(array);
	return(sum);
}

int numColumnsNotNedOutFromTo(struct sequenceMatrix *aSeqMat, int start, int end){
        int i, numSites;
        stringWrap *set, *array;

        set = stringWrapNew(aSeqMat->sampleSize);
        array = stringWrapNew(aSeqMat->sampleSize);
	numSites = 0;

        for( i = start; i < end; i++){
                sequenceMatrix_siteArrayClean(aSeqMat,array,i);
                if (array->length > 1)
                {
                    numSites += 1;
                }
        }

        stringWrapFree(set);
        stringWrapFree(array);

        return numSites;
}

//Fay's Theta_H -- currently broken as NEED TO ADD POLARIZATION
double thetaHPre(stringWrap** set, stringWrap** array, int size)
{
	int i, j;
	double pi, p1, n, nnm1, sum;

	pi = 0.0;
	sum = 0.0;

	for(i=0; i<size; i++)
	{
		if(set[i]->length > 1)
		{
			n = (double)(array[i]->length);
			nnm1 = n * (n-1.0);
			p1 = 0.0;
			pi = 0.0;

			for(j=0; j<set[i]->length; j++)
			{
				p1 = (double)(stringWrapCountChar(array[i], set[i]->cString[j])) / n;
				pi += p1*p1;
			}

			sum += (2.0*pi) / nnm1;
		}
	}

	return sum;
}
//this needs an outgorup to unfold SFS. Assuming that the ancestral sequence is avail
double thetaHFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *ancMat, int start, int end){
	int i;
	double p1, n, nnm1,sum ;
	stringWrap *set, *array, *ancSet;
	char derAllele;

	sum =0.0;
	set = stringWrapNew(aSeqMat->sampleSize);
	ancSet= stringWrapNew(ancMat->sampleSize);
	array = stringWrapNew(aSeqMat->sampleSize);
	for( i = start; i < end; i++){
		sequenceMatrix_siteSetClean(aSeqMat,set,i);
		sequenceMatrix_siteSetClean(ancMat,ancSet,i);
		//two alleles in ingroup and ancRecon at site?
		if(set->length == 2 && ancSet->length == 1){
			sequenceMatrix_siteArrayClean(aSeqMat,array,i);
			//find derivedAllele
			derAllele = (ancSet->cString[0] == set->cString[0]) ? set->cString[1] : set->cString[0];
			n = (double) array->length;
			nnm1 = n * (n - 1.0);
			p1 = (double) stringWrapCountChar(array,derAllele);
			sum += (p1*p1) / nnm1;
		}
	}
	stringWrapFree(set);
	stringWrapFree(ancSet);
	stringWrapFree(array);
	return(2.0*sum);
}

double thetaH(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *ancMat){
	return(thetaHFromTo(aSeqMat,ancMat,0,aSeqMat->length));
}

//Waterson's theta stuff
//
int segSiteCountPre(stringWrap** set, int size)
{
	int i, count;

	count = 0;
	for(i=0; i<size; i++)
	{
		if(set[i]->length > 1)
			count++;
	}

	return count;
}

int segSiteCountFromTo(struct sequenceMatrix *aSeqMat, int start, int end){
	int i,count;
	stringWrap *set;

  	count = 0;
	set = stringWrapNew(aSeqMat->sampleSize);
  	for( i = start; i < end; i++){
		sequenceMatrix_siteSetClean(aSeqMat,set,i);
		if(set->length > 1){
			count+=1;
  		}
	}
	stringWrapFree(set);
	return(count);
}

int segSiteCount(struct sequenceMatrix *aSeqMat){
	return(segSiteCountFromTo(aSeqMat,0,aSeqMat->length));
}

void segSiteLocationsFromTo(struct sequenceMatrix *aSeqMat, vector *segs, int start, int end){
	int i;
	stringWrap *set;
	set = stringWrapNew(aSeqMat->sampleSize + 1);
	vectorInit(segs);
  	for( i = start; i < end; i++){
		sequenceMatrix_siteSetClean(aSeqMat,set,i);
		if(set->length > 1){
			vectorAppend(segs,(double) i);
  		}
	}
	stringWrapFree(set);
}

void segSiteLocations(struct sequenceMatrix *aSeqMat, vector *segs){
	segSiteLocationsFromTo(aSeqMat,segs,0,aSeqMat->length);
}

double thetaWPre(stringWrap** set, int size, int samplesize)
{
	int s;
	s = segSiteCountPre(set, size);
	return (double)(s) / a1f(samplesize);
}

double thetaWFromTo(struct sequenceMatrix *aSeqMat, int start, int end){
	int s;	
	s = segSiteCountFromTo(aSeqMat,start,end);
	return((double)s/a1f(aSeqMat->sampleSize));
}

double thetaW(struct sequenceMatrix *aSeqMat){
	return(thetaWFromTo(aSeqMat,0,aSeqMat->length));
}

//wrapper for tajd function below
double tajD(struct sequenceMatrix *aSeqMat){
	return(tajd(aSeqMat->sampleSize, segSiteCount(aSeqMat), nucdiv(aSeqMat)));
}

double tajDPre(stringWrap** set, stringWrap** array, int size, int samplesize)
{
	double pi = nucdivPre(set, array, size);
	int     s = segSiteCountPre(set, size);
	return tajd(samplesize, s, pi);
}

double tajDFromTo(struct sequenceMatrix *aSeqMat,int start, int end){
	double pi=nucdivFromTo(aSeqMat,start,end);
	int s=segSiteCountFromTo(aSeqMat,start,end);
	return(tajd(aSeqMat->sampleSize,s,pi));
}

double frequency( stringWrap *array, char aChar){
	int i;
	double count=0.0;
  for( i=0; i<array->length; i++) count += ( array->cString[i] == aChar ? 1.0: 0.0 ) ;
  return( count/array->length);
}        

int sampleSizeCleanIndex(struct sequenceMatrix *aSeqMat, int index){
	stringWrap *array;
	int count = 0;
	array = stringWrapNew(aSeqMat->sampleSize);
	sequenceMatrix_siteArrayClean(aSeqMat,array,index);
	count = array->length;
	stringWrapFree(array);
	return(count);
}

double freqAlleleSite(struct sequenceMatrix *aSeqMat, char aChar, int index){
	stringWrap *array;
	double count;
	
	array = stringWrapNew(aSeqMat->sampleSize);
	sequenceMatrix_siteArrayClean(aSeqMat,array,index);
	count = stringWrapCountChar(array,aChar);
	stringWrapFree(array);
	return(count/array->length);
}

//returns most frequency allele
char majorAlleleSite(struct sequenceMatrix *aSeqMat, int index){
	stringWrap *set, *array;
	int i, current , max, tmp;
	char c;
	
	current = max = 0;
	set= stringWrapNew(aSeqMat->sampleSize);
	array = stringWrapNew(aSeqMat->sampleSize);
	sequenceMatrix_siteSetClean(aSeqMat, set,index);
	sequenceMatrix_siteArrayClean(aSeqMat, array,index);
	for(i=0;i<set->length;i++){
		tmp = stringWrapCountChar(array,set->cString[i]);
		if(max < tmp){
			max = tmp;
			current = i;
		}
	}
	c = set->cString[current];
	stringWrapFree(set);
	stringWrapFree(array);
	return(c);
}

//returns most frequency allele
double majorAlleleFreq(struct sequenceMatrix *aSeqMat, int index){
	stringWrap *set, *array;
	int i, current , max, tmp;
	double p;
	
	current = max = 0;
	set= stringWrapNew(aSeqMat->sampleSize);
	array = stringWrapNew(aSeqMat->sampleSize);
	sequenceMatrix_siteSetClean(aSeqMat, set,index);
	sequenceMatrix_siteArrayClean(aSeqMat, array,index);
	for(i=0;i<set->length;i++){
		tmp = stringWrapCountChar(array,set->cString[i]);
		if(max < tmp){
			max = tmp;
			current = i;
		}
	}
	p = (double) max / array->length;
	stringWrapFree(set);
	stringWrapFree(array);
	return(p);
}

//returns least frequency allele
char minorAlleleSite(struct sequenceMatrix *aSeqMat, int index){
	stringWrap *set, *array, *tmpSW;
	int i, current , min, tmp;
	char c;
	
	current = 0;
	min = aSeqMat->sampleSize;
	set= stringWrapNew(aSeqMat->sampleSize);
	array = stringWrapNew(aSeqMat->sampleSize);
	sequenceMatrix_siteSetClean(aSeqMat, set,index);
	sequenceMatrix_siteArrayClean(aSeqMat, array,index);
	for(i=0;i<set->length;i++){
		tmp = stringWrapCountChar(array,set->cString[i]);
		if(min > tmp){
			min = tmp;
			current = i;
		}
	}
	c = set->cString[current];
	//make sure this isn't the same allele if freq = 0.5
	if(((float)min/array->length) == 0.5){
		tmpSW = stringWrapFindOther(set, majorAlleleSite(aSeqMat,index));
		c = tmpSW->cString[0];
		stringWrapFree(tmpSW);
	}
	stringWrapFree(set);
	stringWrapFree(array);
	return(c);
}

int min_rec (struct sequenceMatrix *aSeqMat, int x){
 // Calculate min # rec. events 
  	int a, b, c, e, gtest, flag = 0;
	int size, nsam;
	vector *locs;
	char majAlleleA, majAlleleB;
	
	size = segSiteCount(aSeqMat);
	nsam = aSeqMat->sampleSize;
	locs = vectorNew(size);
	segSiteLocations(aSeqMat,locs);

	if (size<2 || x >= (size-1))
		return (0);
	for (a=x+1; a<size; ++a) {
		for (b=x; b<a; ++b) {
			gtest = 0;
			majAlleleA = majorAlleleSite(aSeqMat,(int)vectorGetIndex(locs,a));
			majAlleleB = majorAlleleSite(aSeqMat,(int)vectorGetIndex(locs,b));
			for (e=0; e<nsam; ++e)
			if (aSeqMat->matrix[e]->cString[b] != majAlleleB && aSeqMat->matrix[e]->cString[a] != majAlleleA) {
				++gtest;
				break;
			}
			for (e=0; e<nsam; ++e)
			if (aSeqMat->matrix[e]->cString[b] != majAlleleB && aSeqMat->matrix[e]->cString[a] == majAlleleA) {
				++gtest;
				break;
			}
			for (e=0; e<nsam; ++e)
			if (aSeqMat->matrix[e]->cString[b] == majAlleleB && aSeqMat->matrix[e]->cString[a] != majAlleleA) {
				++gtest;
				break;
			}
			for (e=0; e<nsam; ++e)
			if (aSeqMat->matrix[e]->cString[b] == majAlleleB && aSeqMat->matrix[e]->cString[a] == majAlleleA) {
				++gtest;
				break;
			}       
			if (gtest == 4) {
				flag = 1;
				break;
			}
		}
		if (flag == 1)
			break;
	}
	if (a==size)
		return (0);
	else {
		c = min_rec(aSeqMat, a);
		return (1+c);
	}
}

double rSquared(struct sequenceMatrix *aSeqMat, int indexA, int indexB){
	char majAlleleA, majAlleleB;
	double pA, pB, pAB, d, denom;
	int comp,i;
	
	majAlleleA = majorAlleleSite(aSeqMat,indexA);
	majAlleleB = majorAlleleSite(aSeqMat,indexB);
//	pA = freqAlleleSite(aSeqMat,majAlleleA,indexA);
//	pB = freqAlleleSite(aSeqMat,majAlleleB,indexB);
	denom = comp = 0;
	pAB =pA = pB = 0.0;
	for(i = 0; i < aSeqMat->sampleSize; i++){
		if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) != 'N' && sequenceMatrix_rowColumn(aSeqMat,i,indexB) != 'N'){
			comp += 1;
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA){
				pA += 1;
			}
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexB) == majAlleleB){
				pB += 1;
			}
			if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA && sequenceMatrix_rowColumn(aSeqMat,i,indexB) == majAlleleB){
				pAB += 1;
			}
		}
	}

	pAB /= comp;
	pA /= comp;
	pB /= comp;
	d = pAB - (pA * pB);
	denom = sqrt(pA * (1.0 - pA) * pB * (1.0-pB));
	d /= denom;
	return(d*d);
}

double rSquaredOmega(struct sequenceMatrix *aSeqMat, int indexA, int indexB){
	char majAlleleA, majAlleleB;
	double pA, pB, pAB, d, denom,sign;
	int comp,i;
	
	majAlleleA = majorAlleleSite(aSeqMat,indexA);
	majAlleleB = majorAlleleSite(aSeqMat,indexB);
//	pA = freqAlleleSite(aSeqMat,majAlleleA,indexA);
//	pB = freqAlleleSite(aSeqMat,majAlleleB,indexB);
	denom = comp = 0;
	pAB =pA = pB = 0.0;
	for(i = 0; i < aSeqMat->sampleSize; i++){
		if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) != 'N' && sequenceMatrix_rowColumn(aSeqMat,i,indexB) != 'N'){
			comp += 1;
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA){
				pA += 1;
			}
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexB) == majAlleleB){
				pB += 1;
			}
			if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA && sequenceMatrix_rowColumn(aSeqMat,i,indexB) == majAlleleB){
				pAB += 1;
			}
		}
	}

	pAB /= comp;
	pA /= comp;
	pB /= comp;
	d = pAB - (pA * pB);
	if(d > 0){
		sign = 1.0;
	}
	else{
		sign = -1.0;
	}
	denom = sqrt(pA * (1.0 - pA) * pB * (1.0-pB));
	d /= denom;
	return (sign * (d*d));
}

void rSquaredCounts(struct sequenceMatrix *aSeqMat, int indexA, int indexB){
	char majAlleleA, majAlleleB;
	double pA, pB, pAB;
	int comp,i;
	
	majAlleleA = majorAlleleSite(aSeqMat,indexA);
	majAlleleB = majorAlleleSite(aSeqMat,indexB);
	comp = 0;
	pAB =pA = pB = 0.0;
	for(i = 0; i < aSeqMat->sampleSize; i++){
		if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) != 'N' && sequenceMatrix_rowColumn(aSeqMat,i,indexB) != 'N'){
			comp += 1;
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA){
				pA += 1;
			}
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexB) == majAlleleB){
				pB += 1;
			}
			if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA && sequenceMatrix_rowColumn(aSeqMat,i,indexB) == majAlleleB){
				pAB += 1;
			}
		}
	}

	pAB /= comp;
	pA /= comp;
	pB /= comp;
	printf("%f\t%f\t%f\t%d",pAB,pA,pB,comp);
}

//jointHeterozygosity returns p(1-p)q(1-q) the squared denom of r^2
double jointHeterozygosity(struct sequenceMatrix *aSeqMat, int indexA, int indexB){
	char majAlleleA, majAlleleB;
	double pA, pB, pAB, d, denom;
	int comp,i;
	
	majAlleleA = majorAlleleSite(aSeqMat,indexA);
	majAlleleB = majorAlleleSite(aSeqMat,indexB);
//	pA = freqAlleleSite(aSeqMat,majAlleleA,indexA);
//	pB = freqAlleleSite(aSeqMat,majAlleleB,indexB);
	denom = comp = 0;
	pAB =pA = pB = 0.0;
	for(i = 0; i < aSeqMat->sampleSize; i++){
		if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) != 'N' && sequenceMatrix_rowColumn(aSeqMat,i,indexB) != 'N'){
			comp += 1;
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA){
				pA += 1;
			}
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexB) == majAlleleB){
				pB += 1;
			}
		}
	}

	pA /= comp;
	pB /= comp;
	d = pAB - (pA * pB);
	denom = pA * (1.0 - pA) * pB * (1.0-pB);
	return(denom);
}

//dij for ZnS
double dij(struct sequenceMatrix *aSeqMat, int indexA, int indexB){
	char majAlleleA, majAlleleB;
	double pA, pB, pAB, d, denom, comp;
	int i;
	
	majAlleleA = majorAlleleSite(aSeqMat,indexA);
	majAlleleB = majorAlleleSite(aSeqMat,indexB);
	pA = freqAlleleSite(aSeqMat,majAlleleA,indexA);
	pB = freqAlleleSite(aSeqMat,majAlleleB,indexB);
	denom = comp = 0;
	pAB =pA = pB = 0.0;
	for(i = 0; i < aSeqMat->sampleSize; i++){
		if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) != 'N' && sequenceMatrix_rowColumn(aSeqMat,i,indexB) != 'N'){
			comp += 1;
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA){
				pA += 1;
			}
			if(sequenceMatrix_rowColumn(aSeqMat,i,indexB) == majAlleleB){
				pB += 1;
			}
			if (sequenceMatrix_rowColumn(aSeqMat,i,indexA) == majAlleleA && sequenceMatrix_rowColumn(aSeqMat,i,indexB) == majAlleleB){
				pAB += 1;
			}
		}
	}
	if(comp < 1){
		return(0);
	}
	else{
		pA /= comp;
		pB /= comp;
		pAB /= comp;
		d = pAB - (pA * pB);
		denom = pA * (1.0 - pA) * pB * (1.0-pB);
		if(denom == 0)
			return(0);
		return((d*d)/denom);
	}
}

double ZnSFromTo(struct sequenceMatrix *aSeqMat, int start, int stop){
	int i, j, size;
	vector *locs;
	double sum= 0.0;
	
	size = segSiteCountFromTo(aSeqMat,start,stop);
	if(size < 2)
		return(0.0);
	locs = vectorNew(size);
	segSiteLocationsFromTo(aSeqMat,locs,start,stop);
	for(i = 0 ; i < size - 1;i++){
		for(j = i+1; j < size; j++){
			sum += dij(aSeqMat,vectorGetIndex(locs,i),vectorGetIndex(locs,j));
			//printf("test: %f\n",sum);
		}
	}
	vectorFree(locs);
	return ((2.0 / (double)(size * (size-1))) * sum);
}

double omegaFromTo(struct sequenceMatrix *aSeqMat, int start, int stop, int site, double **dijTable){
        int i,j, segsites, s;
        vector *locs;
        double sum,sumL,sumR,comp,denom;
        double numer;

        segsites = segSiteCountFromTo(aSeqMat,start,stop);
        s = segsites;
        sum = sumL = sumR =comp=denom=0;
        if(segsites < 3)
                return 0.0;

        locs = vectorNew(segsites);
        segSiteLocationsFromTo(aSeqMat,locs,start,stop);
        //calculate: 
        // sum for denom-- all pairwise r2
        // sumL and sumR for numerator -- blockwise r2
        for(i=0; i<segsites-1; i++){
                for(j=i+1; j<segsites; j++){
                        comp = dijTable[i][j];
                        if(i < site && j >= site)sum += comp;
                        if(i < site && j < site) sumL += comp;
                        if(i >= site && j >= site) sumR += comp;
                }
        }
        denom = sum * (1.0/(site*(s-site)));
        numer = 1.0 / ((site*(site-1)/2) + ((s-site)*(s-site-1)/2));
        numer *= (sumL+sumR);
	vectorFree(locs);
        return(numer/denom);
}

double rEHHatFocalSnp(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *ancMat, int start, int stop, int focalSnpPos)
{
    int i, j, k, l;
    int allsame, currSnpPos, lPos, foundFocalSnp, segsites;
    vector *locs;
    char jBase, kBase, ancBase, derBase;
    int aCount, cCount, gCount, tCount;
    int sameCountDer, sameCountAnc;
    int totalCountDer, totalCountAnc;
    sameCountDer = sameCountAnc = totalCountDer = totalCountAnc = foundFocalSnp = 0;

    assert(ancMat->sampleSize == 1);
    ancBase = sequenceMatrix_rowColumn(ancMat, 0, focalSnpPos);
    if (ancBase == 'A' || ancBase == 'C' || ancBase == 'G' || ancBase == 'T')
    {
        segsites = segSiteCountFromTo(aSeqMat,start,stop);
        locs = vectorNew(1000);
        segSiteLocationsFromTo(aSeqMat, locs, start, stop);
        for(i=0;i<segsites;i++)
        {
            currSnpPos = vectorGetIndex(locs, i);
            if (currSnpPos == focalSnpPos)
                foundFocalSnp = 1;
        }
        //if (!foundFocalSnp){
        //    return NAN;
        //}
        //figure out what the derived allele is
        aCount = cCount = gCount = tCount = 0;
        for(j=0; j < aSeqMat->sampleSize; j++)
        {
            jBase = sequenceMatrix_rowColumn(aSeqMat, j, focalSnpPos);
            if (jBase == 'A' && jBase != ancBase)
                aCount++;
            if (jBase == 'C' && jBase != ancBase)
                cCount++;
            if (jBase == 'G' && jBase != ancBase)
                gCount++;
            if (jBase == 'T' && jBase != ancBase)
                tCount++;
        }
        if (aCount > cCount && aCount > gCount && aCount > tCount)
            derBase = 'A';
        else if (cCount > aCount && cCount > gCount && cCount > tCount)
            derBase = 'C';
        else if (gCount > aCount && gCount > cCount && gCount > tCount)
            derBase = 'G';
        else
            derBase = 'T';
        for(j=0; j < aSeqMat->sampleSize; j++)
        {
            if (sequenceMatrix_rowColumn(aSeqMat, j, focalSnpPos) == derBase)
            {
                for(k=0;k<j;k++)
                {
                    if (sequenceMatrix_rowColumn(aSeqMat, k, focalSnpPos) == derBase)
                    {
                        allsame = 1;
                        //see if j and k are the same at all segsites in the range
                        for(l=0; l<segsites; l++)
                        {
                            lPos = vectorGetIndex(locs, l);
                            if (lPos != focalSnpPos)
                            {
                                jBase = sequenceMatrix_rowColumn(aSeqMat, j, lPos);
                                kBase = sequenceMatrix_rowColumn(aSeqMat, k, lPos);
                                if (jBase != 'N' && kBase != 'N' && jBase != kBase)
                                {
                                    allsame = 0;
                                    break;
                                }
                            }
                        }
                        if (allsame)
                        {
                            sameCountDer++;
                        }
                        totalCountDer++;
                    }
                }
            }
            else if (sequenceMatrix_rowColumn(aSeqMat, j, focalSnpPos) == ancBase)
            {
                for(k=0;k<j;k++)
                {
                    if (sequenceMatrix_rowColumn(aSeqMat, k, focalSnpPos) == ancBase)
                    {
                        allsame = 1;
                        //see if j and k are the same at all segsites in the range
                        for(l=0; l<segsites; l++)
                        {
                            lPos = vectorGetIndex(locs, l);
                            if (lPos != focalSnpPos)
                            {
                                jBase = sequenceMatrix_rowColumn(aSeqMat, j, lPos);
                                kBase = sequenceMatrix_rowColumn(aSeqMat, k, lPos);
                                if (jBase != 'N' && kBase != 'N' && jBase != kBase)
                                {
                                    allsame = 0;
                                    break;
                                }
                            }
                        }
                        if (allsame)
                        {
                            sameCountAnc++;
                        }
                        totalCountAnc++;
                    }
                }
            }
        }
    }
    return (sameCountDer/(float)totalCountDer) / (sameCountAnc/(float)totalCountAnc);
}

/*Get omega at known fixation position.  I.E. center of simulated chromosome*/

double omegaAtCenter(struct sequenceMatrix *aSeqMat, int start, int stop, double site){
	int i,j, segsites, s, loc_i, loc_j;
	vector *locs;
	double sum,sumL,sumR,comp,denom,denomL,denomR;
	double numer;
	//int temp_site = 0;

	segsites = segSiteCountFromTo(aSeqMat,start,stop);
	s = segsites;
	sum = sumL = sumR = comp = denom = denomL = denomR = 0;
	if(segsites < 3)
		return 0.0;
		
	locs = vectorNew(1000);
	segSiteLocationsFromTo(aSeqMat,locs,start,stop);
	
	/*Get snp that is nearest to our fixation*/

	/*min = (fabs(site - vectorGetIndex(locs, 0)));
	for (i = 0; i < segsites; ++i){
		temp_site = vectorGetIndex(locs, i);
		if (fabs(site - temp_site) <= min){
			siteIdx = i;
			min = fabs(site-temp_site);
			}
	}
	*/
	for(i=0; i<segsites-1; i++){
		loc_i = vectorGetIndex(locs, i);
		for(j=i+1; j<segsites; j++){
			loc_j = vectorGetIndex(locs, j);
			comp = dij(aSeqMat, loc_i, loc_j);
			if(loc_i < site && loc_j >= site){
				denom += 1;
				sum += comp;
			}
			if(loc_i < site && loc_j < site){
				denomL += 1;
				sumL += comp;
			}
			if(loc_i >= site && loc_j >= site){
				denomR += 1;
				sumR += comp;		
			}
		}
	}
	denom = sum * (1.0/(double)denom);
	numer = 1.0 / ((double)denomL + (double)denomR);
	numer *= (sumL+sumR);

	if (isnan(denom)){
		return 0.0;
	}
	else{
		return(numer/denom);	
	}
	
}

/*omegaMax -- goes through all possible site divisions to maximize omega
// Kim and Nielsen (2003)
*/

double omegaMaxFromTo(struct sequenceMatrix *aSeqMat, int start, int stop){
        int l, i, j, loc_i, loc_j, segsites;
        double max= 0;
        double tmp=0;
        vector *locs;
        double **dijTable;

	segsites = segSiteCountFromTo(aSeqMat,start,stop);
        dijTable = (double **)malloc( sizeof(double)*segsites );
	locs = vectorNew(segsites);
	segSiteLocationsFromTo(aSeqMat,locs,start,stop);

        if(segsites < 3)
                return(0);
        for(i=0; i<segsites-1; i++){
                dijTable[i] = (double *) malloc( sizeof(double)*segsites );
		loc_i = vectorGetIndex(locs, i);
                for(j=i+1; j<segsites; j++){
                        loc_j = vectorGetIndex(locs, j);
                        dijTable[i][j] = dij(aSeqMat, loc_i, loc_j);
                }
        }
        for(l=3;l<segsites-2;l++){
		tmp = omegaFromTo(aSeqMat,start,stop,l,dijTable);
                if(tmp > max){
                        max = tmp;
                }
        }

	for(i=0; i < segsites-1; i++){
                free(dijTable[i]);
        }
        free(dijTable);
	vectorFree(locs);

        return(max);
}

//outputPolyMask -- writes mask files named for start position
// for use with sweepCoal2 etc.
void outputPolyMask(struct sequenceMatrix *aSeqMat, int start, int stop, FILE *outfile){
	int i, j, seqNumber, nstart, nstop, flag, action;

	nstart = nstop = action = flag = 0;
	seqNumber = 0;
	for(i=0;i<aSeqMat->sampleSize;i++){
		nstart = nstop = action = flag = 0;
		for(j = start; j < stop; j++){
			if (sequenceMatrix_rowColumn(aSeqMat,i,j) == 'N'){
				if (flag){
					nstop++;
				}
				else{
					nstart = j;
					nstop = j;
					flag = 1;
					if (action == 0){
						action = 1;
					}
				}
			}
			else{
				if (flag && action){
					fprintf(outfile,"%d %lf %lf\n",seqNumber,(float)nstart/aSeqMat->length, (float)nstop/aSeqMat->length);
					flag = 0;
				}
			}
		}
		if (flag && action){
			fprintf(outfile,"%d %lf %lf\n",seqNumber,(float)nstart/aSeqMat->length, (float)nstop/aSeqMat->length);
			flag = 0;
		}
		seqNumber++;
	}
	//finish
	fprintf(outfile,"\n//\n");
}

//outputPolyMaskBed -- same as above but uses start and stop to determine normalization
void outputPolyMaskBed(struct sequenceMatrix *aSeqMat, int start, int stop, FILE *outfile){
	int i, j, seqNumber, nstart, nstop, flag, action;

	nstart = nstop = action = flag = 0;
	seqNumber = 0;
	for(i=0;i<aSeqMat->sampleSize;i++){
		nstart = nstop = action = flag = 0;
		for(j = start; j < stop; j++){
			if (sequenceMatrix_rowColumn(aSeqMat,i,j) == 'N'){
				if (flag){
					nstop++;
				}
				else{
					nstart = j;
					nstop = j+1;
					flag = 1;
					if (action == 0){
						action = 1;
					}
				}
			}
			else{
				if (flag && action){
					fprintf(outfile,"%d %lf %lf\n",seqNumber,((float)(nstart-start))/(stop-start), (float)(nstop-start)/(stop-start));
					flag = 0;
				}
			}
		}
		if (flag && action){
			fprintf(outfile,"%d %lf %lf\n",seqNumber,((float)(nstart-start))/(stop-start), (float)(nstop-start)/(stop-start));
			flag = 0;
		}
		seqNumber++;
	}
	//finish
	fprintf(outfile,"\n//\n");
}

/******************** two seqMat things ************/
double fstFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, struct sequenceMatrix *merge, int start,int stop){
	double pi1, pi2, piTot, piW, f;
	pi1 = nucdivFromTo(aSeqMat,start,stop);
	pi2 = nucdivFromTo(bSeqMat,start,stop);
	if(pi1 == 0 && pi2 == 0){
		return(0.0);
	}
	piW = (pi1 * aSeqMat->sampleSize) + (pi2 * bSeqMat->sampleSize);
	piW /= (double) (aSeqMat->sampleSize + bSeqMat->sampleSize);
	piTot = nucdivFromTo(merge,start,stop);
	//printf("test: %f %f %f %f %d\n",pi1,pi2,piW,piTot,merge->sampleSize);
	//exit(1);
	//sequenceMatrix_free(merge);
	f = (piTot - piW) / piTot;
	/*if (f<0)  //Allows the return of negative fsts
		return(0.0);
	else */	
		return(f);
}


//Snn -- Hudson's Snn statistic from 
double SnnFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, int start,int stop){
	double count = 0;
	int i;
	
	for(i=0;i<aSeqMat->sampleSize; i++){
		count += xij_SnnFromTo(aSeqMat,bSeqMat,i,1,start,stop);
	}
	for(i=0;i<bSeqMat->sampleSize; i++){
		count += xij_SnnFromTo(aSeqMat,bSeqMat,i,2,start,stop);
	}
	count /= (double) (aSeqMat->sampleSize + bSeqMat->sampleSize);
	return(count);
}

//xij counts the number of nearest neighbors in same population
double xij_SnnFromTo(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat, \
	int seqIndex1, int belongFlag, int start,int stop){
	
	int i, n1, n2 , nTot;
	double minWith = 667.0;
	double minBet  = 667.0;
	double tmp;
	int minCountW = 0;
	int minCountB = 0;
	
	n1 = aSeqMat->sampleSize;
	n2 = bSeqMat->sampleSize;
	nTot = n1 + n2;
	if(belongFlag == 1){
		//get within min and count
		for(i = 0; i < n1; i++){
			if(i != seqIndex1){
				tmp = seqDist_SnnFromTo(aSeqMat->matrix[seqIndex1],aSeqMat->matrix[i],start,stop);
				if( tmp< minWith){
					minCountW = 1;
					minWith = tmp;
				}
				if( tmp == minWith)
					minCountW += 1;
			}
		}
		//get between min and count
		for(i = 0; i < n2; i++){
			tmp = seqDist_SnnFromTo(aSeqMat->matrix[seqIndex1],bSeqMat->matrix[i],start,stop);
			if( tmp< minBet){
				minCountB = 1;
				minBet = tmp;
			}
			if( tmp == minBet)
				minCountB += 1;
		}		
	}
	else{
		//get within min and count
		for(i = 0; i < n2; i++){
			if(i != seqIndex1){
				tmp = seqDist_SnnFromTo(bSeqMat->matrix[seqIndex1],bSeqMat->matrix[i],start,stop);
				if( tmp< minWith){
					minCountW = 1;
					minWith = tmp;
				}
				if( tmp == minWith)
					minCountW += 1;
			}
		}
		//get between min and count
		for(i = 0; i < n2; i++){
			tmp = seqDist_SnnFromTo(bSeqMat->matrix[seqIndex1],aSeqMat->matrix[i],start,stop);
			if( tmp< minBet){
				minCountB = 1;
				minBet = tmp;
			}
			if( tmp == minBet)
				minCountB += 1;
		}	
	}
	if(minWith < minBet){
		return(1.0);
	}
	if(minWith == minBet){
		return(minCountW / (double) (minCountW+minCountB));
	}
	return(0);

}

double seqDist_SnnFromTo( stringWrap *seq1,  stringWrap *seq2, int start, int stop){
	
	int i;
	double compareCount, nCount;
	double count = 0.0;
	char c1, c2;
	
	compareCount = 0.0;
	nCount = 0.0;
	
	for(i = start; i < stop; i++){
		c1 = seq1->cString[i];
		c2 = seq2->cString[i];
		if(c1 == 'N' || c1 == 'N'){
			nCount += 1; //uncertainty about states?
		}
		else{
			if(c1 != c2)
				count += 1;
		}
		compareCount +=1;
	}
	//arbitrary coverage requirement?
	if (nCount / compareCount > 0.5){
		//return large number
		return(666.0);
	}
	return(count);
}
///////////
////
///
/// still need to convert all below to stringWrap paradigm
///

////////////////////
/*   thetah - pi   */

/*
double hfay( int nsam, int segsites, char **list){
  int s, frequency( char, int, int, char**);
  double pi, p1, nd, nnm1  ;
  
  pi = 0.0 ;
  nd = nsam;
  nnm1 = nd/(nd-1.0) ;
  for( s = 0; s <segsites; s++){
    p1 = frequency('1', s,nsam,list)/nd ;
    pi += 2.0*p1*(2.*p1 - 1.0 )*nnm1 ;
  }
  return( -pi ) ;
}

*/


/*fixed diffs  */
int frequencyFD( char allele,int site,int nsam,  char **list){
  int i, count=0;
  for( i=1; i<nsam; i++){
    count += ( list[i][site] == allele ? 1: 0 );
  }
  return( count);
} 

int fixedDiffs(int segsites, int nsam, char **list){
  int i, fd=0;
  char allele;

  for(i=0; i < segsites; i++){
    allele = list[0][i];
    fd += ((frequencyFD(allele, i, nsam, list) == 0) ? 1:0);
  }
  return(fd);
}

int frequencyFDAllele1( char allele,int site,int nsam,  char **list){
  int count=0;
  count += ( list[1][site] == allele ? 1: 0 ) ;
  return( count);
} 

int fixedDiffsAllele1(int segsites, int nsam, char **list){
  int i, fd=0;
  char allele;

  for(i=0; i < segsites; i++){
    allele = list[0][i];
    fd += (frequencyFDAllele1(allele, i, nsam, list) == 0 ? 1:0);
  }
  return(fd);
}

int ingroupSegSites(int segsites, int nsam, char **list){
  int i, ss = 0;
  
  for(i=0; i < segsites; i++){
    ss += (((frequencyFD('1', i, nsam, list) > 0) && (frequencyFD('1', i, nsam, list) < nsam - 1 )) ? 1:0);
  }
  return(ss);
}

double nucdivIn( int nsam, int segsites, char **list){
  int s;
  double pi, p1, nd, nnm1  ;

  pi = 0.0 ;

  nd = nsam - 1;
  nnm1 = nd/(nd-1.0) ;
  for( s = 0; s <segsites; s++){
    p1 = frequencyFD('1', s,nsam,list)/nd ;
    pi += 2.0*p1*(1.0 -p1)*nnm1 ;
  }
  return( pi ) ;
}

/*haploCount-- returns number of unique haplotypes.
  WARNGING-- sorts seqs in place */
int haploCount(char **list, int nsam){
  int i, count = 1;
  char *tmpSeq;

  
  qsort(list, nsam, sizeof(char *), cmpr); 
  tmpSeq = list[0];
  for(i = 0; i < nsam; i++){
    if (strcmp(list[i],tmpSeq)){
      count += 1;
      tmpSeq = list[i];
    }
  }
  return(count);
}

int cmpr(const void *a, const void *b) { 
 return strcmp(*(char **)a, *(char **)b);
}





	
/************************* tajima.c *************************************************************
 This program calculates Tajima's D when number of sequences, number of segregating sites,
   and average pairwise differences (pi) are known.  It also reports all the coefficients for Tajima's
   D (a1, a2, b1, b2, c1, c2, e1, e2). 
**************************************************************************************************/
double tajd(int nsam, int segsites, double sumk){

  double  a1, a2, b1, b2, c1, c2, e1, e2; 
  
  if( segsites == 0 ) return( 0.0) ;
  
  a1 = a1f(nsam);
  a2 = a2f(nsam);
  b1 = b1f(nsam);
  b2 = b2f(nsam);
  c1 = c1f(a1, b1);
  c2 = c2f(nsam, a1, a2, b2);
  e1 = e1f(a1, c1);
  e2 = e2f(a1, a2, c2);

  return( (sumk - (segsites/a1))/sqrt((e1*segsites) + ((e2*segsites)*(segsites
								      -1))) ) ;
}

double a1f(int nsam){
  double a1;
  int i;
  a1 = 0.0;
  for (i=1; i<=nsam-1; i++) a1 += 1.0/i;
  return (a1);
}


double a2f(int nsam) {
  double a2;
  int i;
  a2 = 0.0;
  for (i=1; i<=nsam-1; i++) a2 += 1.0/(i*i);
  return (a2);
}


double b1f(int nsam){
  double b1;
  b1 = (nsam + 1.0)/(3.0*(nsam-1.0));
  return (b1);
}


double b2f(int nsam){
  double b2;
  b2 = (2*(nsam*nsam + nsam + 3.0))/(9*nsam*(nsam - 1));
  return (b2);
}


double e1f(double a1, double c1){
  double e1;
  e1 = c1/a1;
  return (e1);
}

double e2f(double a1, double a2, double c2){ 
  double e2;
  e2 = c2/((a1*a1)+a2);
  return (e2);
}

double c1f(double a1, double b1){
  double c1;
  c1 = b1 - (1/a1);
  return (c1);
}

double c2f(int nsam, double a1, double a2, double b2){
  double c2;
  c2 = b2 - ((nsam+2)/(a1*nsam)) + (a2/(a1 * a1));
  return (c2);
}


//very general ptr swap
void swap( void **p1,  void **p2)
{
  void *pt = *p1;
  *p1 = *p2;
  *p2 = pt;
}

//for sorting stuff
int compare_doubles(const void *a,const void *b){
	double *pa = (double *) a;
	double *pb = (double *) b;
	if ((*pa - *pb) > 0 ){
		return 1;
	}
	else{
		return - 1;
	}
}

