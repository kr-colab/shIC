import os,sys,random

#usage e.g.: python permuteOrderOfSweepCallRunsByChr.py CEU_100kb.bed 1000 permutedResults/CEU/

realElementsBedFileName, numPermutations, outDir = sys.argv[1:]
numPermutations = int(numPermutations)

def randomlyPickTypeAndPop(runLenLists, disallowed=[]):
    totalLen = 0
    for key in runLenLists:
        if not key in disallowed:
            currLen = len(runLenLists[key])
            totalLen += currLen
    if totalLen == 0:
        return -1, -1

    randDraw = random.random()
    cumProb = 0
    for key in runLenLists:
        if not key in disallowed:
            currLen = len(runLenLists[key])
            prob = currLen / float(totalLen)
            cumProb += prob
            if randDraw <= cumProb:
                return key, runLenLists[key].pop()
    raise Exception # should never get here

def shuffleRunLenLists(runLenLists):
    shuffledRunLensAndTypes = []
    for key in runLenLists:
        random.shuffle(runLenLists[key])
    currType, currItem = randomlyPickTypeAndPop(runLenLists)
    assert not currItem == -1
    shuffledRunLensAndTypes.append((currType, currItem))
    while sum([len(runLenLists[x]) for x in runLenLists]) > 0:
        prevType = currType
        currType, currItem = randomlyPickTypeAndPop(runLenLists, disallowed=prevType)
        if currType == -1:
            assert len(runLenLists[prevType]) > 0
            currType = prevType
            shuffledRunLensAndTypes.append((currType, runLenLists[currType].pop()))
        else:
            shuffledRunLensAndTypes.append((currType, currItem))
    return shuffledRunLensAndTypes

def permuteCoordsFromRunLens(allWinsInChr, runLenLists, c):
    permutedElements = {}
    currRunLenRemaining = 0
    runLenListIndex = 0
    runLenList = shuffleRunLenLists(runLenLists)
    remainingWinsInChr = []
    for s, e, className in allWinsInChr:
        s, e = int(s), int(e)
        if not currRunLenRemaining:
            #initialize the next run
            currRunType, currRunLenRemaining = runLenList[runLenListIndex]
            runLenListIndex += 1
        #continue the current run
        if not currRunType in permutedElements:
            permutedElements[currRunType] = []
        permutedElements[currRunType].append((c, s, e))
        currRunLenRemaining -= 1
    return permutedElements

def appendRunLenForType(runType, runLen, runLensAndTypes):
    if not runType in runLensAndTypes:
        runLensAndTypes[runType] = []
    runLensAndTypes[runType].append(runLen)

def getRunLengthDists(allWins):
    first = True
    runLensAndTypes = {}
    totalRunLen = 0
    for s, e, currType in allWins:
        if first:
            prevType = currType
            runLen = 1
            first = False
        else:
            if prevType == currType:
                runLen += 1
            else:
                appendRunLenForType(prevType, runLen, runLensAndTypes)
                totalRunLen += runLen
                runLen = 1
                prevType = currType
    appendRunLenForType(prevType, runLen, runLensAndTypes)
    totalRunLen += runLen
    assert totalRunLen == len(allWins)
    return runLensAndTypes

def elementNameToType(elementName):
    return elementName.split("_")[0]

def readBedCoordsAsHoLs(bedFileName):
    coordHoLs = {}
    classCountH = {}
    bedFile = open(bedFileName)
    for line in bedFile.xreadlines():
        if not (line.startswith("#") or line.startswith("track")):
            c, s, e, elementName = line.strip().split()[:4]
            elementType = elementNameToType(elementName)
            if not elementType in classCountH:
                classCountH[elementType] = 0
            classCountH[elementType] += 1
            s,e = int(s)+1, int(e)
            if not c in coordHoLs:
                 coordHoLs[c] = []
            coordHoLs[c].append((s, e, elementType))
    bedFile.close()
    return coordHoLs, classCountH

def correctNumPermutedElements(elements, classCountH):
    same = True
    for className in classCountH:
        if classCountH[className] != len(elements[className]):
            same = False

    return same

def writePermutedSet(allElements, classNameLs, outDir, i):
    os.system("mkdir %s/permutation_%s" %(outDir, i))
    for className in classNameLs:
        outFileName = "%s/permutation_%s/%s.bed" %(outDir, i, className)
        with open(outFileName, "w") as outFile:
            for c, s, e in allElements[className]:
                outFile.write("%s\t%s\t%s\n" %(c, s-1, e))

def extendPermutedElements(allPermutedElements, permutedElements):
    for className in permutedElements:
        if not className in allPermutedElements:
            allPermutedElements[className] = []
        allPermutedElements[className].extend(permutedElements[className])

realElements, classCountH = readBedCoordsAsHoLs(realElementsBedFileName)
chrs = list(realElements.keys())
chrs.sort()
for i in range(numPermutations):
    allPermutedElements = {}
    for c in chrs:
        runLensAndTypes = getRunLengthDists(realElements[c])
        permutedElements = permuteCoordsFromRunLens(realElements[c], runLensAndTypes, c)
        extendPermutedElements(allPermutedElements, permutedElements)
    assert correctNumPermutedElements(allPermutedElements, classCountH)
    writePermutedSet(allPermutedElements, sorted(classCountH), outDir, i)
    #sys.stderr.write("done %s permutations-------------\r" %(i+1))
sys.stderr.write("\nall done!\n")
