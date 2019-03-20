import sys, os
from sklearn.externals import joblib
import matplotlib.pyplot as plt
import numpy as np
import random
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn import cross_validation
from scipy.stats import randint as sp_randint
from sklearn import svm
from time import time
from operator import itemgetter
from sklearn.grid_search import GridSearchCV, RandomizedSearchCV
from sklearn import preprocessing

trainingSetDir, classifierPickleFileName = sys.argv[1:3]
statsToUse = sys.argv[3:]

classList = []
trainingData = []
labelToClassName = {}
headerH = {}
for trainingSetFileName in os.listdir(trainingSetDir):
    classList.append(trainingSetFileName.split(".fvec")[0])
    trainingSetFile = open(trainingSetDir + "/" + trainingSetFileName)
    currTrainingData = trainingSetFile.readlines()
    trainingSetFile.close()

    trainingData += currTrainingData[1:]#append all training data from the current set (minus the header)

    currLabelH = {}
    for example in currTrainingData[1:]:
        currLabelH[example.split("\t")[0]] = 1
    assert len(currLabelH) == 1
    labelToClassName[currLabelH.keys()[0]] = trainingSetFileName.split(".fvec")[0]

    header = currTrainingData[0].strip().split("\t")
    headerH[currTrainingData[0].strip()] = 1
    assert header[0] == "classLabel"
    statIndices = []
    if "all" in statsToUse:
        statIndices = range(1, len(header))
    else:
        for i in range(1, len(header)):
            if header[i] in statsToUse or header[i].split("_win")[0] in statsToUse:
                statIndices.append(i)
assert len(headerH) == 1

sys.stderr.write("using these features: %s (indices: %s)\n" %(str(statsToUse), str(statIndices)))
XH = {}
for i in range(len(trainingData)):
    trainingData[i] = trainingData[i].strip().split("\t")
    currVector = []
    if not "nan" in trainingData[i]:
        for j in statIndices:
            currVector.append(float(trainingData[i][j]))
        assert len(currVector) == len(statIndices)
        if not XH.has_key(trainingData[i][0]):
            XH[trainingData[i][0]] = []
        XH[trainingData[i][0]].append(currVector)

#balance the training set
minClassSize = min([len(XH[classLabel]) for classLabel in  XH.keys()])
X = []
y = []
for classLabel in sorted(XH.keys()):
    random.shuffle(XH[classLabel])
    for i in range(minClassSize):
        currVector = XH[classLabel][i]
        X.append(currVector)
        y.append(classLabel)
sys.stderr.write("training set size after balancing: %s\n" %(len(y)))

# Utility function to report best scores
def report(grid_scores, n_top=3):
    top_scores = sorted(grid_scores, key=itemgetter(1), reverse=True)[:n_top]
    for i, score in enumerate(top_scores):
        print("Model with rank: {0}".format(i + 1))
        print("Mean validation score: {0:.3f} (std: {1:.3f})".format(
              score.mean_validation_score,
              np.std(score.cv_validation_scores)))
        print("Parameters: {0}".format(score.parameters))
        print("")


sys.stderr.write("Checking accuracy when distinguishing among all %s classes\n" %(len(XH.keys())))

maxMaxFeatures = len(X[0])
param_grid_forest = {"max_depth": [3, 10, None],
              "max_features": [1, 3, int(maxMaxFeatures**0.5), maxMaxFeatures],
              "min_samples_split": [2, 3, 10],
              "min_samples_leaf": [1, 3, 10],
              "bootstrap": [True, False],
              "criterion": ["gini", "entropy"]}

clf, mlType, paramGrid = ExtraTreesClassifier(n_estimators=100), "extraTreesClassifier", param_grid_forest

heatmap = []
sys.stderr.write("Training %s\n" %(mlType))
grid_search = GridSearchCV(clf,param_grid=param_grid_forest,cv=10,n_jobs=10)
start = time()
grid_search.fit(X, y)
sys.stderr.write("GridSearchCV took %.2f seconds for %d candidate parameter settings.\n"
      % (time() - start, len(grid_search.grid_scores_)))
print "Results for %s" %(mlType)
report(grid_search.grid_scores_)
joblib.dump((labelToClassName, classList, statsToUse, header, grid_search), classifierPickleFileName)
