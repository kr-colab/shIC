CC = gcc
CFLAGS = -O3 -Wall -lm
PG_BASICS =  miscCode/stringWrap.c miscCode/sequenceMatrix.c miscCode/pgSummaryStats.c miscCode/vector.c miscCode/nrutil.c miscCode/bedFile.c

all: niceStats pgStatsBed

niceStats:	niceStats.c miscCode/msGeneralStats.c
	$(CC) niceStats.c miscCode/msGeneralStats.c -o niceStats $(CFLAGS)

pgStatsBed:	pgStatsBedMaskSitesWithMissingData.c $(PG_BASICS)
	$(CC) pgStatsBedMaskSitesWithMissingData.c $(PG_BASICS) -o pgStatsBed $(CFLAGS)

clean:
	$(RM) niceStats pgStatsBed
