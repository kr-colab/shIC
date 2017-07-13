import sys

#python bedSummaryStatsToFeatureVec.py /san/data/1000Genomes/chrX_20kb.ss 5000 0.0 11 0 > /san/data/1000Genomes/chrX_20kb.fvec
#python bedSummaryStatsToFeatureVec.py /san/data/1000Genomes/chr2_20kb_134M_140M_CEU.ss 5000 0.0 11 134000000 > featureVecsToClassify/chr2_20kb_134M_140M_CEU.fvec
bedSummaryStatFileName, minNumSites, piCutoff, numWins, coordOffset = sys.argv[1:]
coordOffset = int(coordOffset) #should be zero unless we are looking at a smaller chunk of a chromosome
minNumSites = int(minNumSites) #going to use 5000 (or 50%) for this for at least the human x chromosome
numWins = int(numWins) #going with 11 here for everything
if numWins % 2 == 0:
    sys.exit("Need to use an odd number of windows!!! AAAARRRRGGGHHHH!!!!\n")
piCutoff = float(piCutoff) #going with 0 (i.e. you get thrown out if you have no SNPs)

def processLineAndCheckPi(line, header, numSites, piCutoff):
    newline = []
    goodLine = False
    thetaH, pi = False, False
    for i in range(len(line)):
        if header[i] == "pi":
            pi = float(line[i])/numSites
            if pi > piCutoff:
                goodLine = True
            newline.append(str(pi))
        elif header[i] == "segSites":
            newline.append(str(float(line[i])/numSites))
        elif header[i] == "thetaH":
            thetaH = float(line[i])/numSites
            newline.append(str(thetaH))
        elif header[i] == "fayWuH":
            if thetaH == False or pi == False:
                raise Exception
            newline.append(str(thetaH-pi))
        else:
            newline.append(line[i])
    if goodLine or pi == False:
        return newline
    else:
        return False

first = True
indicesToKeep = []
linesToKeep = []
bigHeader = ["chrom", "chromStart", "chromEnd"]
bedSummaryStatFile = open(bedSummaryStatFileName)
for line in bedSummaryStatFile.xreadlines():
    #chrom   chromStart      chromEnd        numSites        pi      segSites        thetaH  tajD    fayWuH  HapCount        H1      H12     H2/H1   Omega   ZnS    piWin0  piWin1  piWin2  piWin3  piWin4  piWin5  piWin6  piWin7  piWin8
    if first:
        first = False
        header = line.strip().split("\t")
        for i in range(4, len(header)):
            if not header[i].startswith("piWin"):
                for j in range(numWins):
                    if header[i] == "segSites":
                        bigHeader.append("ss_win%s" %(j))
                    else:
                        bigHeader.append("%s_win%s" %(header[i], j))
                indicesToKeep.append(i)
    else:
        line = line.strip().split("\t")
        c, s, e, numSites = line[:4]
        numSites = int(numSites)
        if numSites >= minNumSites:
            newLine = processLineAndCheckPi(line, header, numSites, piCutoff)
            if newLine:
                linesToKeep.append(newLine)
bedSummaryStatFile.close()
print "\t".join(bigHeader)

def processWin(currWin):
    global numWins, coordOffset
    assert len(currWin) == numWins

    statVecs = {}
    for i in indicesToKeep:
        statVecs[i] = []

    cH = {}
    starts, ends = [], []
    centralC, centralS, centralE = currWin[len(currWin)/2][:3]
    centralS, centralE = int(centralS), int(centralE)
    for subWin in currWin:
        c, s, e, = subWin[:3]
        cH[c] = 1
        starts.append(int(s))
        ends.append(int(e))
        for i in indicesToKeep:
            statVecs[i].append(float(subWin[i]))
    if len(cH.keys()) != 1:
        sys.exit("Coordinates from more than one chromosome found in input file! AAAARRRRGHHHHH!!\n")
    winStart = min(starts)
    winEnd = max(ends)

    outline = [centralC, str(centralS+coordOffset), str(centralE+coordOffset)]
    for i in indicesToKeep:
        minVal = min(statVecs[i])
        if minVal < 0:
            statVecs[i] = [statVal-minVal for statVal in statVecs[i]]
        denom = sum(statVecs[i])
        if denom == 0:
            outline.append("\t".join([str(1.0/len(statVecs[i])) for statVal in statVecs[i]]))
        else:
            outline.append("\t".join([str(statVal/denom) for statVal in statVecs[i]]))
    print "\t".join(outline)

currWin = linesToKeep[:numWins]
for i in range(numWins, len(linesToKeep), 1):
    processWin(currWin)
    currWin.pop(0)
    currWin.append(linesToKeep[i])
processWin(currWin)
