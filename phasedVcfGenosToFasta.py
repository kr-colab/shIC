import sys, os, gzip

try:
    popVcfFileName, refFaFileName,  targetChr, rangeStart, rangeEnd, popFaOutFileName = sys.argv[1:]
    rangeStart, rangeEnd = int(rangeStart), int(rangeEnd)
except Exception:
    usageStr="""phasedVcfGenosToFasta.py -- convert a population sample's phased genotypes in .vcf format to .fasta format
usage:

python popVcfFileName, refFaFileName,  targetChr, rangeStart, rangeEnd, popFaOutFileName

arguments:

popVcfFileName: the path to a VCF file containing phased genotypes for our sample (can be gzipped).
refFaFileName: the path to a fasta file containing the reference sequence for our target chromosome only.
targetChr: the name of the chromosome that we wish to examine (used for labeling fasta entries only).
rangeStart and rangeEnd: the bounds of the region on targetChr that we wish to examine.
popFaOutFileName: path to our output fasta, which includes monomorphic sites and has two entries per diploid.

Note: To mask sites you must edit reference fasta file to have 'N' at sites you wish to ignore.
Another note: This script assumes that all genotypes are imputed (ewww) as in the 1000 Genomes data.
Final note: To use a subset of individuals in your VCF currently you must make a new file omitting the rest.
"""
    sys.exit(usageStr)

sys.stderr.write("reading genotypes\n")

def parseSnpHaps(phasedGenoLs, ref, alt):
    genos = []
    for phasedGeno in phasedGenoLs:
        for allele in phasedGeno.split("|"):
            if allele == '0':
                genos.append(ref)
            else:
                genos.append(alt)
    return genos

if popVcfFileName.endswith(".gz"):
    fopen = gzip.open
else:
    fopen = open

first = True
snpH = {}
lenH = {}
with fopen(popVcfFileName) as popFile:
    for line in popFile:
        if not line.startswith("#"):
            ##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG00096
            line = line.strip().split("\t")
            c, pos, rsid, ref, alt, qual, filt, info, frmt = line[:9]
            pos = int(pos)
            if pos >= rangeStart:
                if pos > rangeEnd:
                    break
                if ref.lower() in ['a','c','g','t'] and alt.lower() in ['a','c','g','t']:
                    snpH[pos] = parseSnpHaps(line[9:], ref, alt)
                    lenH[len(snpH[pos])] = 1
        elif line.startswith("#CHROM"):
            line = line.strip().split("\t")
            sampleIds = line[9:]
assert len(lenH) == 1
numChromos = lenH.keys()[0]
hapIds = []
for sampleId in sampleIds:
    hapIds.append(sampleId+"_1")
    hapIds.append(sampleId+"_2")
assert numChromos == len(hapIds)

refChr = ""
refFile = open(refFaFileName)
sys.stderr.write("reading reference chromosome\n")
for line in refFile.xreadlines():
    if not line.startswith(">"):
        refChr += line.strip().upper()
refFile.close()
toIgnore = {}
for i in range(rangeStart-1, rangeEnd):
    if refChr[i] == 'N':
        toIgnore[i+1] = 1

sys.stderr.write("writing population sample sequence in fasta format\n")
popFaOutFile = open(popFaOutFileName, "w")
for i in range(numChromos):
    popFaOutFile.write(">%s_%s:%s-%s\n" %(hapIds[i], targetChr, rangeStart, rangeEnd))
    doneSites = 0
    for j in range(rangeStart-1, rangeEnd):
        pos = j+1
        if pos <= rangeEnd and pos >= rangeStart:
            if not toIgnore.has_key(pos):
                if snpH.has_key(pos):
                    popFaOutFile.write(snpH[pos][i])
                else:
                    popFaOutFile.write(refChr[j])
            else:
                popFaOutFile.write("N")
        doneSites += 1
        if doneSites % 60 == 0 or pos == rangeEnd or j == len(refChr)-1:
            popFaOutFile.write("\n")
popFaOutFile.close()
