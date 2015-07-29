/*msGeneralStats.h -- contains general analysis of Hudson format output */

#ifndef MSGS_INC
#define MSGS_INC
#include "stdio.h"


int **imatrix(int nsam,int len);

//code from achaz
double * compute_HarmonicSums( int n );
double beta_FU1995( int i, double *HarmonicSums, int n  );
double sigma_ii_FU1995( int i, double *HarmonicSums, int n );
double sigma_ij_FU1995( int i, int j, double *HarmonicSums, int n );

int frequency( char allele,int site,int nsam,  char **list);
double nucdiv(int, int, char **);
double thetah(int, int, char **);
int segSites(int, int, char **);
void getSiteFreqSpec(int segsites, int nsam, char**list, int nSites, int *derived_counts);
int sortcmp (int *n1, int *n2);
int getHaplotypeFreqSpec(int segsites, int nsam, char **list, int *haplotype_counts);
double meanEHH(int segsites, double *posit, double delta, int nsam, char **list);
double meanREHH(int segsites, double *posit, double delta, int nsam, char **list);
void petrovHStatsWindow(int segsites, int nwins, double *posit, double *winsH12, double *winsH1, double *winsH2, int nsam, char **list);
double petrovH1(int *haplotype_counts, int nsam);
double petrovH2(int *haplotype_counts, int nsam);
double petrovH12(int *haplotype_counts, int nsam);
int nHaplotypes(int segsites, int nsam, char **list);
int biggerlist(int nsam,unsigned nmax,char **list );
int sampleSizeSite(int site, int nsam, char **list);
double ZnS(int segsites, int nsam, char** list);
double dij(int i, int j, int nsam, char** list);
double omega(int left, int segsites,int nsam, char** list);
double omegaWithTable(int left, int segsites, int nsam, char** list, double** dijTable);
double omegaMax(int segsites,int nsam, char** list);
double omegaCenter(int siteIdx , int segsites,int nsam, char** list);
double sigmaIJ(int i, int j, int n);
double sigmaII(int i, int n);
double sigmaBeta(int i,int n);
double sigmaAlpha(int n);
double achazNeutTestExponentWeights(int nsam, int segsites, char **list, int exponent1, int exponent2, double *harmonicSums);
double achazThetaExponentWeights(int nsam, int segsites, char **list, int exponent);
double achazThetaHPi(int nsam, int segsites, char **list);
double achazTajimasDExtreme(int nsam, int segsites, char **list);
void nucdivWindow( int nwins, double *posit, double *output, int nsam, int segsites, char **list);
void tajdWindow(int nwins, double *posit, double *output, int nsam, int segsites, char **list);
void getSiteFreqSpecWindow(int segsites, int nsam, char**list, int nSites, int *derived_counts, double *pos, double low, double high);


//subpopn declares
int frequencySub(char allele, int site, int startAllele, int stopAllele, char **list);
double nucdivSub( int nsam, int segsites, int startAllele, int stopAllele, char **list);
double thetahSub( int nsam, int segsites, int startAllele, int stopAllele, char **list);
int nHaplotypesSub(int segsites, int nsam, int startAllele, int stopAllele, char **list);
double ZnSSub(int segsites, int nsam, int startAllele, int stopAllele, char** list);
double fst2Subs(int segsites, int nsam, int start1, int stop1, int start2, int stop2, char **list);
int segSitesSub(int segsites, int nsam, int startAllele, int stopAllele, char **list);
double Snn(int segsites,int nsam, int n1, int n2, char **list);
double xij_Snn(int segsites,int nsam, int seqIndex1, int n1, int n2, char **list);
double seqDist_Snn(int segsites, int index1,  int index2, char **list);
double meanRefDist(int segsites, int nsam, char **list);

void nucdivSubWindow( int nwins, double *posit, double *output, int nsam, int segsites,int startAllele, int stopAllele, char **list);
void fst2SubsWindow(int nwins, double *posit, double *output,int segsites, int nsam, int start1, int stop1, int start2, int stop2, char **list);


//twoSite
void sampleConfig(int i, int j, int nsam, char** list, int *config);
void printPairwiseSampleConfigs(int segSites, int nsam, char **list,double *posit, int nsites);
void printPairwiseSampleConfigs2Popn(int segsites, int nsam, int popnSize1, char **list, double *posit, int nsites);
void sampleConfig2Popn(int i, int j, int nsam, int popnSize1, char** list, int *config);

//From Hudson
double tajd(int nsam, int segsites, double sumk);
double a1f(int);
double a2f(int);
double b1f(int);
double b2f(int);
double c1f(double, double);
double c2f(int, double, double, double);
double e1f(double, double);
double e2f(double, double, double);

#endif
