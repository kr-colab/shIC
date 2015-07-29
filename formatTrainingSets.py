import sys, os, random

neutTrainingFileName, softTrainingFilePrefix, hardTrainingFilePrefix, sweepWin, outDir = sys.argv[1:]
sweepWin = int(sweepWin)

sweepFilePaths, linkedFilePaths = {}, {}
for trainingFilePrefix in [softTrainingFilePrefix, hardTrainingFilePrefix]:
    trainingSetDir = "/".join(trainingFilePrefix.split("/")[:-1])
    trainingFilePrefixDirless = trainingFilePrefix.split("/")[-1]
    linkedFilePaths[trainingFilePrefix] = []
    sweepFilePaths[trainingFilePrefix] = []

    for fileName in os.listdir(trainingSetDir):
        if fileName.startswith(trainingFilePrefixDirless):
            #equilibSoft_1.msout
            winNum = int(fileName.split("_")[1].split(".")[0])
            if winNum == sweepWin:
                sweepFilePaths[trainingFilePrefix].append(trainingSetDir + "/" + fileName)
            else:
                linkedFilePaths[trainingFilePrefix].append(trainingSetDir + "/" + fileName)

def getExamplesFromFVFile(simFileName):
    simFile = open(simFileName)
    lines = [line.strip() for line in simFile.readlines()]
    header = lines[0]
    examples = lines[1:]
    simFile.close()
    return header, examples

def getExamplesFromFVFileLs(simFileLs):
    examples = []
    for filePath in simFileLs:
        header, currExamples = getExamplesFromFVFile(filePath)
        examples += currExamples
    return examples

header, neutExamples = getExamplesFromFVFile(neutTrainingFileName)
linkedSoftExamples = getExamplesFromFVFileLs(linkedFilePaths[softTrainingFilePrefix])
softExamples = getExamplesFromFVFileLs(sweepFilePaths[softTrainingFilePrefix])
linkedHardExamples = getExamplesFromFVFileLs(linkedFilePaths[hardTrainingFilePrefix])
hardExamples = getExamplesFromFVFileLs(sweepFilePaths[hardTrainingFilePrefix])
trainingSetLs = [linkedSoftExamples, softExamples, linkedHardExamples, hardExamples]
for i in range(len(trainingSetLs)):
    random.shuffle(trainingSetLs[i])
    trainingSetLs[i] = trainingSetLs[i][:len(neutExamples)]
linkedSoftExamples, softExamples, linkedHardExamples, hardExamples = trainingSetLs

outFileNames = ["neut.fvec", "linkedSoft.fvec", "soft.fvec", "linkedHard.fvec", "hard.fvec"]
outExamples = [neutExamples, linkedSoftExamples, softExamples, linkedHardExamples, hardExamples]
for i in range(len(outFileNames)):
    outFile = open(outDir + "/" + outFileNames[i], "w")
    outFile.write("classLabel\t%s\n" %(header))
    for example in outExamples[i]:
        outFile.write("%s\t%s\n" %(i, example))
    outFile.close()
