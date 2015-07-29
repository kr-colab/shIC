import sys, os
import numpy as np
from sklearn.externals import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn import cross_validation
from scipy.stats import randint as sp_randint
from sklearn import svm
from time import time
from operator import itemgetter
from sklearn.grid_search import GridSearchCV, RandomizedSearchCV
from sklearn import preprocessing

#bedGraph example: python classifyWinsInChromosome.py classifiers/shIT/tennessenHumanTest/tennessen_humanTest_ETc.p featureVecsToClassify/chr18_200kb_CEU.fvec 200 bedGraph /home/dan/public_html/spatialSVM/browserStuff/shIT_chr18_200kb_CEU_test

classifierPickleFileName, genomicFeatureVecFileName, winSizeInKb, bedOrBedGraph = sys.argv[1:5]
if not bedOrBedGraph in ["bed","bedGraph"]:
    sys.exit("bedOrBedGraph must be \"bed\" or \"bedGraph\"\n")
elif bedOrBedGraph == "bedGraph":
    bedGraphOutPrefix = sys.argv[5]
    sys.stderr.write("Running classifier stored in %s on chromosome represented in %s; then writing posterior probability estimates in .bedGraph format\n" %(classifierPickleFileName, genomicFeatureVecFileName))
else:
    sys.stderr.write("Running classifier stored in %s on chromosome represented in %s; then printing classifications in .bed format\n" %(classifierPickleFileName, genomicFeatureVecFileName))


labelToClassName, classList, statsToUse, header, grid_search = joblib.load(classifierPickleFileName)

statIndices = []
for i in range(1, len(header)):
    if "all" in statsToUse or header[i] in statsToUse or header[i].split("_win")[0] in statsToUse:
        statIndices.append(i)

genomicFeatureVecFile = open(genomicFeatureVecFileName)
genomeData = genomicFeatureVecFile.readlines()
genomicFeatureVecFile.close()

genomeHeader = genomeData[0].strip().split("\t")
assert genomeHeader[3:] == header[1:]
genomeData = genomeData[1:]#remove the header

newGenomeData = []
for example in genomeData:
    if not "nan" in example:
        newGenomeData.append(example)
genomeData = newGenomeData

genomeX = []
coords = []
for i in range(len(genomeData)):
    genomeData[i] = genomeData[i].strip().split("\t")
    currVector = []
    coords.append(genomeData[i][:3])
    for j in range(len(genomeData[i])):
        if j-2 in statIndices:
            assert genomeHeader[j] == header[j-2]
            currVector.append(float(genomeData[i][j]))
    assert len(currVector) == len(statIndices)
    genomeX.append(currVector)
genomeX = np.array(genomeX)

if bedOrBedGraph == "bed":
    predictions = grid_search.predict(genomeX)
    assert len(coords) == len(predictions)
    sys.stderr.write("Made predictions for %s instances in %s\n" %(len(genomeX), genomicFeatureVecFileName))
    predCounts = {}
    for className in classList:
        predCounts[className] = 0

    classToColorStr = {"hard": "255,0,0", "linkedHard": "255,150,150", "soft": "0,0,255", "linkedSoft": "150,150,255", "neut": "0,0,0"}
    outlines = ["track name=shIClassifications%s description=\"Classification from soft/hard inference tool (%s kb wins)\" visibility=2 itemRgb=On" %(winSizeInKb, winSizeInKb)]
else:
    probas = grid_search.predict_proba(genomeX)
    sys.stderr.write("Estimated posterior probabilities for %s instances in %s\n" %(len(genomeX), genomicFeatureVecFileName))
    predCounts = {}
    assert len(coords) == len(probas)
    outFiles = []
    classNameToDisplayName = {"hard": "Hard", "linkedHard": "Hard-linked", "soft": "Soft", "linkedSoft": "Soft-linked", "neut": "Neutral"}
    for classLabel in range(len(labelToClassName.keys())):
        className = labelToClassName[str(classLabel)]
        outFiles.append(open("%s.%s.probs.bedGraph" %(bedGraphOutPrefix, className), "w"))
        outFiles[-1].write("track type=bedGraph name=shIClassificationProbs%s%s description=\"Posterior for %s class from soft/hard inference tool (%s kb wins)\" visibility=full color=0,150,0 autoScale=off viewLimits=0.0:1.0\n" %(classNameToDisplayName[className], winSizeInKb, classNameToDisplayName[className], winSizeInKb))
    
for i in range(len(coords)):
    c, s, e = coords[i]
    if bedOrBedGraph == "bed":
        predictedClass = labelToClassName[predictions[i]]
        predCounts[predictedClass] += 1
        outlines.append("%s\t%s\t%s\t%s_%s_%s\t0\t.\t%s\t%s\t%s" %(c, s, e, predictedClass, c, predCounts[predictedClass], s, e, classToColorStr[predictedClass]))
    else:
        for j in range(len(probas[i])):
            outFiles[j].write("%s\t%s\t%s\t%s\n" %(c, s, e, probas[i][j]))

if bedOrBedGraph == "bed":
    for classLabel in sorted(predCounts):
        sys.stderr.write("predicted %s elements in class %s (%s of total)\n" %(predCounts[classLabel], classLabel, predCounts[classLabel]/float(len(predictions))))

    for outline in outlines:
        print outline
else:
    for i in range(len(outFiles)):
        outFiles[i].close()
