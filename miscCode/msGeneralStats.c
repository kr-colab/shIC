#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "msGeneralStats.h"
/* allocates space for gametes (character strings) */


char **
	cmatrix(nsam,len)
	int nsam, len;
{
	int i;
	char **m;

	if( ! ( m = (char **) malloc( (unsigned)( nsam*sizeof( char* )) ) ) )
		perror("alloc error in cmatrix") ;
	for( i=0; i<nsam; i++) {
		if( ! ( m[i] = (char *) malloc( (unsigned) (len*sizeof( char )) )))
			perror("alloc error in cmatric. 2");
	}
	return( m );
}

int **
	imatrix(nsam,len)
	int nsam, len;
{
	int i,ii;
	int **m;

	if( ! ( m = (int **) malloc( (unsigned)( nsam*sizeof( int* )) ) ) )
		perror("alloc error in imatrix") ;
	for( i=0; i<nsam; i++) {
		if( ! ( m[i] = (int *) malloc( (unsigned) (len*sizeof( int )) )))
			perror("alloc error in imatric. 2");
		for(ii=0;ii<len;ii++)m[i][ii]=0;
	}
	return( m );
}
int
	biggerlist(nsam, nmax, list )
	int nsam ;
unsigned nmax ;
char ** list ;
{
	int i;

	for( i=0; i<nsam; i++){
		list[i] = (char *)realloc( list[i],nmax*sizeof(char) ) ;
		if( list[i] == NULL ) perror( "realloc error. bigger");
	}
	return(0);
}                        


//below we've corrected the sampleSizes for missing data
double nucdiv( int nsam, int segsites, char **list)
{
	int s, frequency( char, int, int, char**);
	double pi, p1, nd, nnm1  ;

	pi = 0.0 ;


	for( s = 0; s <segsites; s++){
		nd = sampleSizeSite(s,nsam,list);
		if (nd > 1){
			nnm1 = nd/(nd-1.0) ;
			p1 = frequency('1', s,nsam,list)/nd ;
			pi += 2.0*p1*(1.0 -p1)*nnm1 ;
		}
	}
	return( pi ) ;
}

//fills a vector size l with values of nucdiv in "windows"
void nucdivWindow( int nwins, double *posit, double *output, int nsam, int segsites, char **list)
{
	int s, frequency( char, int, int, char**);
	int wcount = 0;
	double pi, p1, nd, nnm1  ;
	double start, end, delta;
	start = 0;
	delta = 1.0 / nwins;
	end = delta;
	
	while(start < 1.0){
		pi = 0.0 ;
		for( s = 0; s <segsites; s++){
			if(posit[s] <=end && posit[s] > start){
				nd = sampleSizeSite(s,nsam,list);
				if (nd > 1){
					nnm1 = nd/(nd-1.0) ;
					p1 = frequency('1', s,nsam,list)/nd ;
					pi += 2.0*p1*(1.0 -p1)*nnm1 ;
				}
			}
		}
		output[wcount++]=pi;
		start += delta;
		end += delta;
	}
}


//fills a vector size l with values of Tajima's D in "windows"
void tajdWindow(int nwins, double *posit, double *output, int nsam, int segsites, char **list)
{
        int s, frequency( char, int, int, char**);
        int wcount = 0;
        double pi, p1, nd, nnm1  ;
        double start, end, delta;
        int segsitesInWin;
        start = 0;
        delta = 1.0 / nwins;
        end = delta;

        while(start < 1.0){
                pi = 0.0 ;
                segsitesInWin = 0;
                for( s = 0; s <segsites; s++){
                        if(posit[s] <=end && posit[s] > start){
                                nd = sampleSizeSite(s,nsam,list);
                                if (nd > 1){
                                        segsitesInWin += 1;
                                        nnm1 = nd/(nd-1.0) ;
                                        p1 = frequency('1', s,nsam,list)/nd ;
                                        pi += 2.0*p1*(1.0 -p1)*nnm1 ;
                                }
                        }
                }
                output[wcount++]=tajd(nsam,segsitesInWin,pi);
                start += delta;
                end += delta;
        }
}


double achazThetaExponentWeights(int nsam, int segsites, char **list, int exponent)
{
        int s, frequency( char, int, int, char**);
        double thetaA, i, wi, nd, nnm1, wSum ;

        wSum = 0.0 ;
        for ( i = 1; i < nsam; i++){
            wi = pow(i,exponent);
            wSum += wi;
        }

        thetaA = 0.0 ;
        for( s = 0; s <segsites; s++){
                nd = sampleSizeSite(s,nsam,list);
                if (nd > 1){
                        nnm1 = nd/(nd-1.0);
                        i = frequency('1', s,nsam,list);
                        wi = pow(i,exponent);
                        thetaA += wi*i;
                }
        }

        thetaA = thetaA/wSum;
        return thetaA;
}

//pi with a twist of theta H
double achazThetaHPi(int nsam, int segsites, char **list)
{
        int s, frequency( char, int, int, char**);
        double thetaA, i, wi, nd, nnm1, wSum ;

        wSum = 0.0 ;
        for ( i = 1; i < nsam; i++){
            wi = i*(nsam-i);
            wSum += wi;
        }

        thetaA = 0.0 ;
        for( s = 0; s <segsites; s++){
                nd = sampleSizeSite(s,nsam,list);
                if (nd > 1){
                        nnm1 = nd/(nd-1.0);
                        i = frequency('1', s,nsam,list);
                        wi = i*(nsam-i);
                        thetaA += wi*i;
                }
        }

        thetaA = thetaA/wSum;
        return thetaA;
}

//pi minus a sumary stat that is highest when your SFS is U-shaped
double achazTajimasDExtreme(int nsam, int segsites, char **list)
{
        int s, frequency( char, int, int, char**);
        double pi, upsideDownPi, i, wi1, wi2, nd, nnm1, wSum1, wSum2 ;

        wSum1 = 0.0 ;
        wSum2 = 0.0;
        for ( i = 1; i < nsam; i++){
            wi1 = (nsam-i);
            wi2 = pow(((nsam/2.0)-i),2)/i;
            wSum1 += wi1;
            wSum2 += wi2;
        }

        pi = 0.0 ;
        upsideDownPi = 0.0;
        for( s = 0; s <segsites; s++){
                nd = sampleSizeSite(s,nsam,list);
                if (nd > 1){
                        nnm1 = nd/(nd-1.0);
                        i = frequency('1', s,nsam,list);
                        wi1 = (nsam-i);
                        wi2 = pow(((nsam/2.0)-i),2)/i;
                        pi += wi1*i;
                        upsideDownPi += wi2*i;
                }
        }

        pi = pi/wSum1;
        upsideDownPi = upsideDownPi/wSum2;
        return pi-upsideDownPi;
}

double sigmaAlpha(int n)
{
    int i;
    int sum = 0;
    for (i=1; i<n;i++)
    {
        sum += 1.0/(double)i;
    }
    return sum;
}

double sigmaBeta(int i,int n)
{
    double part1,part2,part3;

    part1 = (double)(2*n)/(double)((n-i+1)*(n-i));
    part2 = sigmaAlpha(n+1) - sigmaAlpha(i);
    part3 = 2.0/(double)(n-i);

    return part1*part2-part3;
}

double sigmaII(int i, int n)
{
    if (2*i < n)
    {
        return sigmaBeta(i+1,n);
    }
    else if (2*i == n)
    {
        return 2.0*(sigmaAlpha(n) - sigmaAlpha(i))/(double)(n-i) - (1.0/(double)(i*i));
    }
    else
    {
        return sigmaBeta(i,n) - (1.0/(double)(i*i));
    }
}

//copied from Achaz
double * compute_HarmonicSums( int n ){
	double *HarmonicSums;
	int i;
	HarmonicSums = (double *)malloc( sizeof(double)*n );
	if(!HarmonicSums)
		fprintf(stderr, "compute_HarmonicSums: malloc error for HarmonicSums\n"), exit(3);
	HarmonicSums[0]=0;
	i=1;
	while(i<n){
		HarmonicSums[i] = HarmonicSums[i-1] + 1.0/i;
		i++;
	}
	return HarmonicSums;
}

//these next 3 functions are copied from Achaz because my versions have a bug somewhere and I am lazy; these seem to work
double beta_FU1995( int i, double *HarmonicSums, int n  ){

	double ai = HarmonicSums[i-1], an = HarmonicSums[n-1];
	double beta=0;
	double nloci= (double)n;

	beta = 2.0 * nloci * ( an + (1.0/nloci) - ai )/(  (nloci-i+1.0 )*(nloci-i) ) - 2.0/(nloci - i);

	return beta;
}

double sigma_ii_FU1995( int i, double *HarmonicSums, int n ){
	double nloci= (double)n;
	double sigma_ii=0;
	double ai = HarmonicSums[i-1], an = HarmonicSums[n-1];
	if( 2*i < n )
	{
		sigma_ii = beta_FU1995( i+1, HarmonicSums, n  );
	}
	else
	{
		if( 2*i == n  )
		{
			sigma_ii = 2.0*(an-ai)/(nloci - i) - 1.0/(i*i);
		}
		else
		{
			sigma_ii = beta_FU1995( i , HarmonicSums, n  ) - 1.0/(i*i);
		}
	}
	return sigma_ii;
}

double sigma_ij_FU1995( int i, int j, double *HarmonicSums, int n ){
	double nloci= (double)n;
	double sigma_ij=0;

	if(i==j){
		return sigma_ii_FU1995( i, HarmonicSums, n );
	}
 	if(i<j){
		int tmp=i;
		i=j;
		j=tmp;
	}
	double  ai=HarmonicSums[i-1],
	       aj=HarmonicSums[j-1],
	       an=HarmonicSums[n-1];
	if( i+j < n )
	{
		sigma_ij = ( beta_FU1995( i+1, HarmonicSums, n  ) -  beta_FU1995( i, HarmonicSums, n  ) ) / 2.0;
	}
	else
	{
		if( i+j == n  )
		{
			sigma_ij  =  ((an - ai)/(nloci - i) + (an - aj)/(nloci - j)) 
			           - ( ( beta_FU1995( i, HarmonicSums, n   ) +  beta_FU1995( j+1 , HarmonicSums, n ) )  / 2.0 ) 
				   - (1.0/(i*j));

		}
		else
		{
			sigma_ij = (( beta_FU1995( j, HarmonicSums, n  ) -  beta_FU1995( j+1, HarmonicSums, n ) )/2.0) - (1.0/(i*j));
		}

	}
	return sigma_ij;

}

double sigmaIJ(int i, int j, int n)
{
    double part1,part2,part3,part4;
    int tmp;

    if (i<j)
    {
        tmp = i;
        i = j;
        j = tmp;
    }
    if (i+j < n)
    {
        return ((sigmaBeta(i+1,n) - sigmaBeta(i,n))/2.0);
    }
    else if (i+j == n)
    {
        part1 = (sigmaAlpha(n) - sigmaAlpha(i)) / (double)(n-i);
        part2 = (sigmaAlpha(n) - sigmaAlpha(j)) / (double)(n-j);
        part3 = (sigmaBeta(i,n) + sigmaBeta(j+1,n)) / 2.0;
        part4 = 1.0/(double)(i*j);
        return part1 + part2 - part3 - part4;
    }
    else
    {
        //printf("%d,%d,sigmaIJ=%f\n",i,j,((sigmaBeta(j,n) - sigmaBeta(j+1,n))/2.0 - 1.0/(double)(i*j)));
        return ((sigmaBeta(j,n) - sigmaBeta(j+1,n))/2.0 - 1.0/(double)(i*j));
    }
}

double achazNeutTestExponentWeights(int nsam, int segsites, char **list, int exponent1, int exponent2, double *harmonicSums)
{
        int i, j, s, frequency( char, int, int, char**);
        double nd, T, Tnum, alphaN, w1sum, w2sum, betaNpart1, betaNpart2, betaN, thetaEstimate, thetaSquaredEstimate, a2;
        int *afs;
        double *w1;
        double *w2;
        double *omega;

        //initializing afs, and weight arrays, and weight sums
        //in each array, element index zero remains uninitialized: not counting monomorphic sites
        afs = (int *) malloc(sizeof(int)*nsam);
        w1 = (double *) malloc(sizeof(double)*nsam);
        w2 = (double *) malloc(sizeof(double)*nsam);
        w1sum = 0.0;
        w2sum = 0.0;
        for (i = 1; i < nsam; i++){
            afs[i] = 0;
            w1[i] = pow(i,exponent1);
            w2[i] = pow(i,exponent2);
            w1sum += w1[i];
            w2sum += w2[i];
        }

        omega = (double *) malloc(sizeof(double)*nsam);
        //initializing omega, again not using the element at index 0
        for (i = 1; i < nsam; i++){
            omega[i] = (w1[i]/w1sum) - (w2[i]/w2sum);
        }

        //compute afs
        for( s = 0; s <segsites; s++){
                nd = sampleSizeSite(s,nsam,list);
                if (nd > 1){
                        i = frequency('1', s,nsam,list);
                        afs[i] += 1;
                }
        }

        //calculate neutrality test statistic
        Tnum = 0.0 ;
        betaNpart1 = 0.0;
        betaNpart2 = 0.0;
        alphaN = 0.0;
        for ( i = 1; i < nsam; i++){

            Tnum += omega[i]*(double)(i*afs[i]);

            alphaN += i*omega[i]*omega[i];
            betaNpart1 += (double)(i*i)*omega[i]*omega[i]*sigma_ii_FU1995(i,harmonicSums,nsam);
            for (j=i+1; j<nsam; j++)
            {
                betaNpart2 += (double)(2*i*j)*omega[i]*omega[j]*sigma_ij_FU1995(i,j,harmonicSums,nsam);
                //printf("omega[i]=%f,omega[j]=%f,betaNpart2: i=%d,j=%dn=%d,sigmaIJ=%f\n", omega[i],omega[j],i,j,nsam,sigma_ij_FU1995(i,j,harmonicSums,nsam));
            }
        }

        //borrowing some code here from Achaz to get estimates of theta (using thetaW) and theta^2
        a2 = 0;
        for(i=1;  i< nsam;  i++){
		a2 += 1.0/(double)(i*i);
	}
	thetaEstimate =  (double) segsites / harmonicSums[nsam-1];
	thetaSquaredEstimate =  (double) segsites*( (double)segsites-1.0 )/(harmonicSums[nsam-1]*harmonicSums[nsam-1] + a2 );

        betaN=betaNpart1+betaNpart2;
        //printf("betaNpart1=%f,betaNpart2=%f\tbetaN=%f\n",betaNpart1,betaNpart2,betaN);
        //printf("alphaN=%f\n",alphaN);
        //printf("root=%f\n",pow(alphaN*thetaEstimate+(betaN*thetaSquaredEstimate),0.5));
        T = Tnum / sqrt(alphaN*thetaEstimate+(betaN*thetaEstimate*thetaSquaredEstimate));

        free(afs);
        free(w1);
        free(w2);
        free(omega);
        return T;
}


/* Fay's theta_H  */
//again corrected for sampleSize variation
double	thetah( int nsam, int segsites, char **list)
{
	int s, frequency( char, int, int, char**);
	double pi, p1, nd  ;

	pi = 0.0 ;

	nd = nsam;
	for( s = 0; s <segsites; s++){
		p1 = frequency('1', s,nsam,list) ;
		nd = sampleSizeSite(s,nsam,list);
		if(nd > 1){
			pi += (p1*p1)/( nd*(nd-1.0) ) ; 
		}
	}
	return(pi*2.0) ;
}


int frequency( char allele,int site,int nsam,  char **list){
	int i, count=0;
	for( i=0; i<nsam; i++) count += ( list[i][site] == allele ? 1: 0 ) ;
	return( count);
}        

//sampleSizeSite -- returns the sampleSize at a site corrected for missing data
int sampleSizeSite(int site, int nsam, char **list){
	return(nsam - frequency('N',site,nsam,list));	
}


//gets the right number of segSites when there are Ns
int segSites(int segsites, int nsam, char **list){
	int i, ss = 0;

	for(i=0; i < segsites; i++){
		ss += (((frequency('1', i, nsam, list) > 0) && (frequency('0', i, nsam, list) >0 )) ? 1:0);
	}
	return(ss);
}

//gets the derived site frequency spectrum; fixations treated as monomorphic
//derived_counts must be an int array of length nsam
//derived_counts[i] is the fraction of sites with derived allele present in i chromosomes
//number of monomorphic sites is sorted in derived_counts[0]
void getSiteFreqSpec(int segsites, int nsam, char**list, int nSites, int *derived_counts)
{
        int i;
        int freq;

        for (i=0; i<nsam; i++)
        {
            derived_counts[i] = 0;
        }

        int polycount = 0;
        for (i=0; i<segsites; i++)
        {
                freq = frequency('1', i, nsam, list);
                if (freq > 0 && freq < nsam)
                {
                        polycount++;
                        derived_counts[freq] += 1;
                }
        }
        derived_counts[0] = nSites-polycount;
}

//gets the derived site frequency spectrum; fixations treated as monomorphic
//derived_counts must be an int array of length nsam
//derived_counts[i] is the fraction of sites with derived allele present in i chromosomes
//number of monomorphic sites is sorted in derived_counts[0]
void getSiteFreqSpecWindow(int segsites, int nsam, char**list, int nSites, int *derived_counts, double *pos, double low, double high)
{
        int i;
        int freq;

        for (i=0; i<nsam; i++)
        {
            derived_counts[i] = 0;
        }

        int polycount = 0;
        for (i=0; i<segsites; i++)
        {
		if (pos[i] >low && pos[i] <= high)
		{
                	freq = frequency('1', i, nsam, list);
                	if (freq > 0 && freq < nsam)
                	{
                        	polycount++;
                        	derived_counts[freq] += 1;
                	}
		}
        }
        derived_counts[0] = nSites-polycount;
}

//counts the number of haplotypes, and gets their frequencies (stored in haplotype_counts)
//haplotype_counts must be an int array of length nsam
//haplotype_counts[i] is the number of haplotypes found in exactly i+1 individuals
int getHaplotypeFreqSpec(int segsites, int nsam, char **list, int *haplotype_counts)
{
        int i;
        int j;
        int k;
        int haplotype_found;
        int allsame;
        int freq;

        int n_haplotypes = 0;
        char haplotypes[nsam][segsites+1];
        int haplotype_occurrences[nsam];

        for(i=0; i<nsam; i++)
        {
            haplotype_counts[i] = 0;
            haplotype_occurrences[i] = 0;
        }

        for(i=0; i<nsam; i++)
        {
                haplotype_found = 0;
                for(j=0; j<n_haplotypes; j++)
                {
                        allsame = 1;
                        for(k=0; k<segsites; k++)
                        {
                                if(haplotypes[j][k] != list[i][k])
                                {
                                        if(haplotypes[j][k] == 'N' )
                                                haplotypes[j][k] = list[i][k];

                                        else if(list[i][k] != 'N')
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
                        for(j=0; j<segsites; j++)
                        {
                                haplotypes[n_haplotypes-1][j] = list[i][j];
                        }
                        haplotypes[n_haplotypes-1][segsites]='\0';
                        haplotype_occurrences[n_haplotypes-1]=1;
                }
        }

        for (i=0; i<n_haplotypes; i++)
        {
                freq = haplotype_occurrences[i];
                if (freq > 0 && freq <= nsam)
                {
                        haplotype_counts[freq-1] += 1;
                }
        }

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

//gets H12, H1, and H2 in windows
void petrovHStatsWindow(int segsites, int nwins, double *posit, double *winsH12, double *winsH1, double *winsH2, int nsam, char **list)
{
	int i;
	int j;
	int k;
	int haplotype_found, freq;
	int allsame;
	float start, delta, end;
	int n_haplotypes = 0;
	char haplotypes[nsam][segsites+1];
        int haplotype_occurrences[nsam];
        int haplotype_counts[nsam];
        int wcount = 0;
	start = 0;
        delta = 1.0 / nwins;
        end = delta;

        while(start < 1.0)
        {
                for(i=0; i<nsam; i++)
                {
                    haplotype_counts[i] = 0;
                    haplotype_occurrences[i] = 0;
                }
                for(i=0; i<nsam; i++)
                {
                    haplotype_found = 0;
                    for(j=0; j<n_haplotypes; j++)
                    {
                        allsame = 1;
                        for(k=0; k<segsites; k++)
                        {
                            if (posit[k] > start && posit[k] <= end)
                            {
                                if(haplotypes[j][k] != list[i][k])
                                {
                                    if(haplotypes[j][k] == 'N' )
                                        haplotypes[j][k] = list[i][k];
                                    else if(list[i][k] != 'N')
                                    {
                                        allsame = 0;
                                        break;
                                    }
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
                        for(j=0; j<segsites; j++)
                        {
                            haplotypes[n_haplotypes-1][j] = list[i][j];
                        }
                        haplotypes[n_haplotypes-1][segsites]='\0';
                        haplotype_occurrences[n_haplotypes-1]=1;
                    }
                }
                for (i=0; i<n_haplotypes; i++)
                {
                    freq = haplotype_occurrences[i];
                    if (freq > 0 && freq <= nsam)
                    {
                        haplotype_counts[freq-1] += 1;
                    }
                }
                winsH12[wcount]=petrovH12(haplotype_counts,nsam);
                winsH1[wcount]=petrovH1(haplotype_counts,nsam);
                winsH2[wcount]=petrovH2(haplotype_counts,nsam);
                wcount++;
                start += delta;
                end += delta;
                n_haplotypes = 0;
        }
}

double meanEHH(int segsites, double *posit, double delta, int nsam, char **list)
{
	int i, j, k, l;
	int allsame;
	float start, end;
        int sameCount;
        int totalCount;
        float ehhSum = 0;
        float ehhCount = 0;
        for(i=0;i<segsites;i++)
        {
            sameCount = 0;
            totalCount = 0;
            start = posit[i]-delta;
            end = posit[i]+delta;
            if (start >= 0 && end <= 1)
            {
                for(j=0; j<nsam; j++)
                {
                    if (list[j][i] == '1')
                    {
                        for(k=0;k<j;k++)
                        {
                            if (list[k][i] == '1')
                            {
                                allsame = 1;
                                //see if j and k are the same at all segsites in the range
                                for(l=0; l<segsites; l++)
                                {
                                    if (posit[l] >= start && posit[l] <= end && l != i)
                                    {
                                        if (list[j][l] != list[k][l])
                                        {
                                            allsame = 0;
                                            break;
                                        }
                                    }
                                }
                                if (allsame)
                                {
                                    sameCount++;
                                }
                                totalCount++;
                            }
                        }
                    }
                }
                if (totalCount > 0)
                {
                    ehhSum += sameCount/(float)totalCount;
                }
                ehhCount += 1;
            }
        }
        return ehhSum/ehhCount;
}

double meanREHH(int segsites, double *posit, double delta, int nsam, char **list)
{
	int i, j, k, l;
	int allsame;
	float start, end;
        int sameCountDer, sameCountAnc;
        int totalCountDer, totalCountAnc;
        float ehhSum = 0;
        float ehhCount = 0;
        for(i=0;i<segsites;i++)
        {
            sameCountDer = 0;
            totalCountDer = 0;
            sameCountAnc = 0;
            totalCountAnc = 0;
            start = posit[i]-delta;
            end = posit[i]+delta;
            if (start >= 0 && end <= 1)
            {
                for(j=0; j<nsam; j++)
                {
                    if (list[j][i] == '1')
                    {
                        for(k=0;k<j;k++)
                        {
                            if (list[k][i] == '1')
                            {
                                allsame = 1;
                                //see if j and k are the same at all segsites in the range
                                for(l=0; l<segsites; l++)
                                {
                                    if (posit[l] >= start && posit[l] <= end && l != i)
                                    {
                                        if (list[j][l] != list[k][l])
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
                    else if (list[j][i] == '0')
                    {
                        for(k=0;k<j;k++)
                        {
                            if (list[k][i] == '0')
                            {
                                allsame = 1;
                                //see if j and k are the same at all segsites in the range
                                for(l=0; l<segsites; l++)
                                {
                                    if (posit[l] >= start && posit[l] <= end && l != i)
                                    {
                                        if (list[j][l] != list[k][l])
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
                if (totalCountAnc > 0 && totalCountDer > 0 && sameCountAnc > 0)
                {
                    ehhSum += (sameCountDer/(float)totalCountDer) / (sameCountAnc/(float)totalCountAnc);
                }
                ehhCount += 1;
            }
        }
        return ehhSum/ehhCount;
}

//counts number of haplotypes
int nHaplotypes(int segsites, int nsam, char **list)
{
	int i;
	int j;
	int k;
	int haplotype_found;
	int allsame;
	
	int n_haplotypes = 0;
	char haplotypes[nsam][segsites+1];
	
	for(i=0; i<nsam; i++)
	{
		haplotype_found = 0;
		for(j=0; j<n_haplotypes; j++)
		{
			allsame = 1;
			for(k=0; k<segsites; k++)
			{
				if(haplotypes[j][k] != list[i][k])
				{
					if(haplotypes[j][k] == 'N' )
						haplotypes[j][k] = list[i][k];
												
					else if(list[i][k] != 'N')
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
			for(j=0; j<segsites; j++)
					haplotypes[n_haplotypes-1][j] = list[i][j];
                        haplotypes[n_haplotypes-1][segsites]='\0';
		}
	}
		
	return n_haplotypes;
}

double dij(int i, int j, int nsam, char** list){
	
	double pi  = 0.0;
	double pj  = 0.0;
	double pij = 0.0;
	double count = 0.0;
	int k;
	for(k=0; k<nsam; k++){
		if(list[k][i] != 'N' && list[k][j] != 'N'){
			if(list[k][i] == '1')
				pi++;
			
			if(list[k][j] == '1')
				pj++;
			
			if(list[k][i] == '1' && list[k][j] == '1')
				pij++;	
			count += 1;
		}
	}
	if (count == 0){
		return(0);
	}
	else{
		pi  /= count;
		pj  /= count;
		pij /= count;
	
		double Dij = pij - (pi*pj);
	
		return (Dij*Dij) / ((pi*(1.0-pi)) * (pj*(1.0-pj)));	
	}
}


/* implements the ZnS statistic of J Kelly (1997), "A Test of Neutrality Based */
/* on Interlocus Associations." Genetics, 146: 1197-1206.                      */
double ZnS(int segsites, int nsam, char** list){

	if(segsites < 2)
		return 0.0;

	int i,j,k, s;
	double sum = 0.;
	
	s = segSites(segsites, nsam, list);
	for(i=0; i<segsites-1; i++){
                if (frequency('1', i, nsam, list) < nsam){
			for(j=i+1; j<segsites; j++){
        	        	if (frequency('1', j, nsam, list) < nsam){
					sum += dij(i, j, nsam, list);
				}
			}
		}
	}
	return (2.0 / (double)(s * (s-1))) * sum;
}

/*omega statistic from Kim and Nielsen (2003)
** not robust to missing data
*/
double omegaWithTable(int left, int segsites,int nsam, char** list, double** dijTable){
	int i,j, s;
	double sum,sumL,sumR,comp,denom;
	double numer;
	
	sum = sumL = sumR =comp=denom=0;
	if(segsites < 3)
		return 0.0;
	s = segSites(segsites, nsam, list);
	//calculate: 
	// sum for denom-- all pairwise r2
	// sumL and sumR for numerator -- blockwise r2
	for(i=0; i<segsites-1; i++){
		for(j=i+1; j<segsites; j++){
			comp = dijTable[i][j];
			if(i < left && j >= left)sum += comp;
			if(i < left && j < left) sumL += comp;
			if(i >= left && j >= left) sumR += comp;	
		}
	}
	denom = sum * (1.0/(left*(s-left)));
	numer = 1.0 / ((left*(left-1)/2) + ((s-left)*(s-left-1)/2));
	numer *= (sumL+sumR);
	//printf("n/d: %f d: %f n %f sumL: %f sumR: %f term: %f left: %d\n",numer/denom,denom,numer,sumL,sumR,1.0 / ((left*(left-1)/2) + ((s-left)*(s-left-1)/2)),left);
	return(numer/denom);
}

/*omega statistic from Kim and Nielsen (2003)
** not robust to missing data
*/
double omega(int left, int segsites,int nsam, char** list){
        int i,j, s;
        double sum,sumL,sumR,comp,denom;
        double numer;

        sum = sumL = sumR =comp=denom=0;
        if(segsites < 3)
                return 0.0;
        s = segSites(segsites, nsam, list);
        //calculate: 
        // sum for denom-- all pairwise r2
        // sumL and sumR for numerator -- blockwise r2
        for(i=0; i<segsites-1; i++){
                for(j=i+1; j<segsites; j++){
                        comp = dij(i, j, nsam, list);
                        if(i < left && j >= left)sum += comp;
                        if(i < left && j < left) sumL += comp;
                        if(i >= left && j >= left) sumR += comp;
                }
        }
        denom = sum * (1.0/(left*(s-left)));
        numer = 1.0 / ((left*(left-1)/2) + ((s-left)*(s-left-1)/2));
        numer *= (sumL+sumR);
        return(numer/denom);
}

/*Get omega at known fixation position. I.E. center of simulation*/
double omegaCenter(int siteIdx , int segsites,int nsam, char** list){
	int i,j, s;
	double sum,sumL,sumR,comp,denom;
	double numer;
	
	sum = sumL = sumR =comp=denom=0;
	if(segsites < 3)
		return 0.0;
	s = segSites(segsites, nsam, list);
	//calculate: 
	// sum for denom-- all pairwise r2
	// sumL and sumR for numerator -- blockwise r2
	for(i=0; i<segsites-1; i++){
		for(j=i+1; j<segsites; j++){
			comp = dij(i, j, nsam, list);
			if(i < siteIdx && j >= siteIdx)sum += comp;
			if(i < siteIdx && j < siteIdx) sumL += comp;
			if(i >= siteIdx && j >= siteIdx) sumR += comp;	
		}
	}
	denom = sum * (1.0/(siteIdx*(s-siteIdx)));
	numer = 1.0 / ((siteIdx*(siteIdx-1)/2) + ((s-siteIdx)*(s-siteIdx-1)/2));
	numer *= (sumL+sumR);
//	printf("d: %f n %f sumL: %f sumR: %f term: %f left: %d\n",denom,numer,sumL,sumR,1.0 / ((left*(left-1)/2) + ((s-left)*(s-left-1)/2)),left);
	if (isnan(denom)){
        return 0.0;
    }
    else{
        return(numer/denom);    
    }
}


/*omegaMax -- goes through all possible site divisions to maximize omega
// Kim and Nielsen (2003)
// This version builds a table of all pairwise r^2 vals for faster downstream computation
*/
double omegaMax(int segsites,int nsam, char** list){
	int l, i, j;
	double max= 0;
	double tmp=0;
	double **dijTable;

	dijTable = (double **)malloc( sizeof(double)*segsites );

	if(segsites < 3)
		return(0);
	for(i=0; i<segsites-1; i++){
		dijTable[i] = (double *) malloc( sizeof(double)*segsites );
		for(j=i+1; j<segsites; j++){
			dijTable[i][j] = dij(i, j, nsam, list);
                        //printf("dijTable[i][j]: %f\n", dijTable[i][j]);
		}
	}
	for(l=3;l<segsites-2;l++){
		tmp = omegaWithTable(l, segsites, nsam,list, dijTable);
		if(tmp > max){
			max = tmp;
		}
	}

	for(i=0; i < segsites-1; i++){
                free(dijTable[i]);
        }
        free(dijTable);

	return(max);
}

/////////////////////////
//Two Site Utils
//////
////
//
//sampleConfig-- fills vector with the sample configuration
// for sites i and j
void sampleConfig(int i, int j, int nsam, char** list, int *config){
	int p1, p2, x11, k;
	p1 = p2 = x11 = 0;

	for(k=0; k<nsam; k++){
		if(list[k][i] != 'N' && list[k][j] != 'N'){
			if(list[k][i] == '1')
				p1++;
			
			if(list[k][j] == '1')
				p2++;
			
			if(list[k][i] == '1' && list[k][j] == '1')
				x11++;	
		}
	}
	if(p1 > p2){
		config[0] = p1;
		config[1] = p2;
	}
	else{
		config[0] = p2;
		config[1] = p1;
	}
	config[2] = x11;
}

void printPairwiseSampleConfigs(int segsites, int nsam, char **list, double *posit, int nsites){
	int i,j, config[3];

	if(segsites < 2){
		return;
	}
	else{
		for(i=0; i<segsites-1; i++){
			for(j=i+1; j<segsites; j++){
				sampleConfig(i, j, nsam, list,config);
				printf("%d %d %d %d %d\n",nsam,config[0],config[1],config[2],(int)floor((posit[j]-posit[i]) * nsites));
			}	
		}
	}
}

//sampleConfig2Popn- fills vector with the sample configuration
// for sites i and j; 6 dimensions for 2 population 2 site sample config
void sampleConfig2Popn(int i, int j, int nsam, int popnSize1, char** list, int *config){
	int p1, p2, x11,p3,p4,y11, k;
	p1 = p2 = x11 = p3 = p4 = y11 = 0;

	for(k=0; k<nsam; k++){
		if(list[k][i] != 'N' && list[k][j] != 'N'){
			if(k < popnSize1){
				if(list[k][i] == '1')
					p1++;
			
				if(list[k][j] == '1')
					p2++;
			
				if(list[k][i] == '1' && list[k][j] == '1')
					x11++;	
			}
			else{
				if(list[k][i] == '1')
					p3++;
			
				if(list[k][j] == '1')
					p4++;
			
				if(list[k][i] == '1' && list[k][j] == '1')
					y11++;
			}
		}
	}
//	if(p1 > p2){
		config[0] = p1;
		config[1] = p2;
//	}
//	else{
//		config[0] = p2;
//		config[1] = p1;
//	}
	config[2] = x11;
//	if(p3 > p4){
		config[3] = p3;
		config[4] = p4;
//	}
//	else{
//		config[3] = p4;
//		config[4] = p3;
//	}
	config[5] = y11;
}

void printPairwiseSampleConfigs2Popn(int segsites, int nsam, int popnSize1, char **list, double *posit, int nsites){
	int i,j, config[6];

	if(segsites < 2){
		return;
	}
	else{
		for(i=0; i<segsites-1; i++){
			for(j=i+1; j<segsites; j++){
				sampleConfig2Popn(i, j, nsam, popnSize1,list,config);
				printf("%d %d %d %d %d %d %d %d %d\n",popnSize1,nsam-popnSize1,config[0],config[1],config[2],config[3],config[4],config[5],(int)floor((posit[j]-posit[i]) * nsites));
			}	
		}
	}
}
/****************
///   Sub popn versions
******************/

//frequencySub-- allows for arbitrary allele indexes as you would need for sub pops
int frequencySub(char allele, int site, int startAllele, int stopAllele, char **list){
	int i, count=0;
	for( i=startAllele; i<stopAllele; i++) count += ( list[i][site] == allele ? 1: 0 ) ;
	return( count);
}


//gets the right number of segSites when there are Ns
int segSitesSub(int segsites, int nsam, int startAllele, int stopAllele, char **list){
	int i, ss = 0;

	for(i=0; i < segsites; i++){
		ss += (((frequencySub('1', i, startAllele, stopAllele, list) > 0) && (frequencySub('0', i, startAllele, stopAllele, list) >0 )) ? 1:0);
	}
	return(ss);
}


//sampleSizeSiteSub -- returns the sampleSize at a site corrected for missing data
//in startAllele to stopAllele rows
int sampleSizeSiteSub(int site, int nsam, int startAllele, int stopAllele, char **list){
	return(stopAllele - startAllele  - frequencySub('N',site,startAllele, stopAllele,list));	
}


double nucdivSub( int nsam, int segsites, int startAllele, int stopAllele, char **list){
	int s;
	double pi, p1, nd, nnm1  ;

	pi = 0.0 ;


	for( s = 0; s <segsites; s++){
		nd = sampleSizeSiteSub(s,nsam,startAllele,stopAllele,list);
		if (nd > 1){
			nnm1 = nd/(nd-1.0) ;
			p1 = frequencySub('1', s,startAllele,stopAllele,list)/nd ;
			pi += 2.0*p1*(1.0 -p1)*nnm1 ;
		}
	}
	return( pi ) ;
}


//fills a vector size l with values of nucdiv in "windows". for subpops
void nucdivSubWindow( int nwins, double *posit, double *output, int nsam, int segsites,int startAllele, int stopAllele, char **list)
{
	int s, frequency( char, int, int, char**);
	int wcount = 0;
	double pi, p1, nd, nnm1  ;
	double start, end, delta;
	start = 0;
	delta = 1.0 / nwins;
	end = delta;
	
	while(start < 1.0){
		pi = 0.0 ;
		for( s = 0; s <segsites; s++){
			if(posit[s] <=end && posit[s] > start){
				nd = sampleSizeSiteSub(s,nsam,startAllele,stopAllele,list);
				if (nd > 1){
					nnm1 = nd/(nd-1.0) ;
					p1 = frequencySub('1', s,startAllele,stopAllele,list)/nd ;
					pi += 2.0*p1*(1.0 -p1)*nnm1 ;
				}
			}
		}
		output[wcount++]=pi;
		start += delta;
		end += delta;
	}
}


void fst2SubsWindow(int nwins, double *posit, double *output,int segsites, int nsam, int start1, int stop1, int start2, int stop2, char **list){
	double h1[nwins], h2[nwins], hTot[nwins], hW;

	int i;
	
	nucdivSubWindow(nwins,posit,h1,nsam, segsites, start1, stop1, list);
	nucdivSubWindow(nwins,posit,h2,nsam, segsites,start2,stop2,list);
	nucdivWindow(nwins,posit,hTot,nsam, segsites,list);
	
	for(i=0;i<nwins;i++){
		hW = (((stop1 - start1) * h1[i]) + ((stop2-start2) * h2[i])) / (double) ((stop1-start1) +(stop2-start2));

		if(h1[i] == 0.0 && h2[i] == 0.0)output[i]=0.0;
		else output[i] = (hTot[i]-hW) / hTot[i];
	}
	
}

/* Fay's theta_H  */
//again corrected for sampleSize variation
double thetahSub( int nsam, int segsites, int startAllele, int stopAllele, char **list){
	int s;
	double pi, p1, nd  ;

	pi = 0.0 ;

	nd = nsam;
	for( s = 0; s <segsites; s++){
		p1 = frequencySub('1', s,startAllele,stopAllele,list) ;
		nd = sampleSizeSiteSub(s,nsam,startAllele,stopAllele,list);
		if(nd > 1 && p1 != (stopAllele-startAllele)){
			pi += (p1*p1)/( nd*(nd-1.0) ) ; 
		}
	}
	return(pi*2.0) ;
}

//counts number of haplotypes
int nHaplotypesSub(int segsites, int nsam, int startAllele, int stopAllele, char **list)
{
	int i;
	int j;
	int k;
	int haplotype_found;
	int allsame;
	
	int n_haplotypes = 0;
	char haplotypes[ stopAllele - startAllele][segsites];

	for(i=startAllele; i<stopAllele; i++)
	{
		haplotype_found = 0;
		for(j=0; j<n_haplotypes; j++)
		{
			allsame = 1;
			for(k=0; k<segsites; k++)
			{
				if(haplotypes[j][k] != list[i][k])
				{
					if(haplotypes[j][k] == 'N' )
						haplotypes[j][k] = list[i][k];
												
					else if(list[i][k] != 'N')
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
			for(j=0; j<segsites; j++)
					haplotypes[n_haplotypes-1][j] = list[i][j];
		}
	}
		
	return n_haplotypes;
}

double dijSub(int i, int j, int nsam, int startAllele, int stopAllele, char** list){
	
	double pi  = 0.0;
	double pj  = 0.0;
	double pij = 0.0;
	double count = 0.0;
	double Dij = 0.0;
	double denom =0.0;
	int k;
	for(k=startAllele; k<stopAllele; k++){
		if(list[k][i] != 'N' && list[k][j] != 'N'){
			if(list[k][i] == '1')
				pi++;
			
			if(list[k][j] == '1')
				pj++;
			
			if(list[k][i] == '1' && list[k][j] == '1')
				pij++;	
			count += 1;
		}
	}
	if (count == 0){
		return(0);
	}
	else{
		pi  /= count;
		pj  /= count;
		pij /= count;
	
		Dij = pij - (pi*pj);
	    denom = (pi*(1.0-pi)) * (pj*(1.0-pj));
		if(denom == 0)
			return(0);
		return ((Dij*Dij) / denom);	
	}
}

/* implements the ZnS statistic of J Kelly (1997), "A Test of Neutrality Based */
/* on Interlocus Associations." Genetics, 146: 1197-1206.                      */
double ZnSSub(int segsites, int nsam, int startAllele, int stopAllele, char** list){

	if(segSitesSub(segsites,nsam,startAllele,stopAllele,list) < 2){
		return(0);
	}

	int i,j, s;
	double sum = 0.;
	
	s = segSitesSub(segsites, nsam, startAllele,stopAllele,list);
	for(i=0; i<segsites-1; i++){
		for(j=i+1; j<segsites; j++){
			sum += dijSub(i, j, nsam,startAllele,stopAllele, list);
		}
	}
	
	return (2.0 / (double)(s * (s-1))) * sum;
}

//fst--
double fst2Subs(int segsites, int nsam, int start1, int stop1, int start2, int stop2, char **list){
	double h1, h2, hTot, hW;
	double f;
	
	h1 = nucdivSub(nsam, segsites, start1, stop1, list);
	h2 = nucdivSub(nsam, segsites,start2,stop2,list);
	hTot = nucdiv(nsam, segsites,list);
	hW = (((stop1 - start1) * h1) + ((stop2-start2) * h2)) / (double) ((stop1-start1) +(stop2-start2));
	//printf("%f %f %f %f\n",h1,h2,hTot,hW);
	if(h1 == 0.0 && h2 == 0.0)
		return(0.0);
	f = (hTot-hW) / hTot;
//	if(f < 0)
//		return(0.0); going to allow this to return negative values... might be better for estimation?
	return(f);
	
}
//******************
//Snn -- Hudson's Snn statistic from Hudson (2000)
double Snn(int segsites,int nsam, int n1, int n2, char **list){
	double count = 0;
	int i;
	
	for(i=0;i<nsam; i++){
		count += xij_Snn(segsites,nsam,i,n1,n2,list);
	}
	count /= (double) nsam;
	return(count);
}

//this counts the actually proportions of nearest neighbors in same pop
double xij_Snn(int segsites,int nsam, int seqIndex1, int n1, int n2, char **list){
	
	int i;
	double minWith = 667.0;
	double minBet  = 667.0;
	double tmp;
	int minCountW = 0;
	int minCountB = 0;
	
	if(seqIndex1 < n1){
		//get within min and count
		for(i = 0; i < n1; i++){
			if(i != seqIndex1){
				tmp = seqDist_Snn(segsites,seqIndex1,i,list);
				if( tmp< minWith){
					minCountW = 1;
					minWith = tmp;
				}
				if( tmp == minWith)
					minCountW += 1;
			}
		}
		//get between min and count
		for(i = n1; i < nsam; i++){
			tmp = seqDist_Snn(segsites,seqIndex1,i,list);
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
		for(i = n1; i < nsam; i++){
			if(i != seqIndex1){
				tmp = seqDist_Snn(segsites,seqIndex1,i,list);
				if( tmp< minWith){
					minCountW = 1;
					minWith = tmp;
				}
				if( tmp == minWith)
					minCountW += 1;
			}
		}
		//get between min and count
		for(i = 0; i < n1; i++){
			tmp = seqDist_Snn(segsites,seqIndex1,i,list);
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

//seqDist_Snn is the metric for Snn
double seqDist_Snn(int segsites, int index1,  int index2, char **list){
	
	int i;
	double compareCount, nCount;
	double count = 0.0;
	char c1, c2;
	
	compareCount = 0.0;
	nCount = 0.0;
	
	for(i = 0; i < segsites; i++){
		c1 = list[index1][i];
		c2 = list[index2][i];
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

//meanRefDist-- calculates the mean dist of all alleles to a ref--sequence 1
double meanRefDist(int segsites, int nsam, char **list){
	int i;
	double sum = 0.0;
	
	for(i = 1; i < nsam; i++){
		sum += seqDist_Snn(segsites,0,i,list);
	}
	return(sum/(nsam-1));
}

//From Hudson's ms package
/************************* tajima.c *************************************************************
 This program calculates Tajima's D when number of sequences, number of segregating sites,
   and average pairwise differences (pi) are known.  It also reports all the coefficients for Tajima's
   D (a1, a2, b1, b2, c1, c2, e1, e2). 
**************************************************************************************************/


	double
tajd(int nsam, int segsites, double sumk)
{

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

double a1f(int nsam)
{
double a1;
int i;
a1 = 0.0;
for (i=1; i<=nsam-1; i++) a1 += 1.0/i;
return (a1);
}


double a2f(int nsam) 
{
double a2;
int i;
a2 = 0.0;
for (i=1; i<=nsam-1; i++) a2 += 1.0/(i*i);
return (a2);
}


double b1f(int nsam)
{
double b1;
b1 = (nsam + 1.0)/(3.0*(nsam-1.0));
return (b1);
}


double b2f(int nsam)
{
double b2;
b2 = (2*(nsam*nsam + nsam + 3.0))/(9*nsam*(nsam - 1));
return (b2);
}


double e1f(double a1, double c1) 
{
double e1;
e1 = c1/a1;
return (e1);
}

double e2f(double a1, double a2, double c2)
{ 
double e2;
e2 = c2/((a1*a1)+a2);
return (e2);
}


double c1f(double a1, double b1)
{
double c1;
c1 = b1 - (1/a1);
return (c1);
}


double c2f(int nsam, double a1, double a2, double b2)
{
double c2;
c2 = b2 - ((nsam+2)/(a1*nsam)) + (a2/(a1 * a1));
return (c2);
}
