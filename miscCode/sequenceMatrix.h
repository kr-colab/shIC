/*sequenceMatrix.h -- simple object to hold sequences */
/* ADK 10/9/09 */

#ifndef SEQMAT_INC
#define SEQMAT_INC

#include "stringWrap.h"
#include <stdlib.h>
#include <stdio.h>
//#include "pgSummaryStats.h"



struct sequenceMatrix{
	int sampleSize, length;
	stringWrap **matrix;
	stringWrap **names;
}sequenceMatrix;

struct sequenceMatrix *sequenceMatrix_new(int sampleSize, int length);
void sequenceMatrix_free(struct sequenceMatrix *aSeqMat);
struct sequenceMatrix *sequenceMatrix_importFasta(char *fileName);

void sequenceMatrix_siteArray(struct sequenceMatrix *aSeqMat, stringWrap *dst, int index);
void sequenceMatrix_siteArrayClean(struct sequenceMatrix *aSeqMat, stringWrap *dst, int index);
void sequenceMatrix_siteSet(struct sequenceMatrix *aSeqMat, stringWrap *dst, int index);
void sequenceMatrix_siteSetClean(struct sequenceMatrix *aSeqMat, stringWrap *dst, int index);
char sequenceMatrix_rowColumn(struct sequenceMatrix *aSeqMat, int r, int c);

void sequenceMatrix_rowColumn_set(struct sequenceMatrix *aSeqMat,int r, int c, char x);

void sequenceMatrix_NOutSitesWithTooMuchMissingData(struct sequenceMatrix *aSeqMat, double maxFracMissingData);
void sequenceMatrix_printFasta(struct sequenceMatrix *aSeqMat);
struct sequenceMatrix *sequenceMatrix_merge(struct sequenceMatrix *aSeqMat, struct sequenceMatrix *bSeqMat);
#endif

