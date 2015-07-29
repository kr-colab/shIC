import sys, os

statFilePrefix = sys.argv[1]

summaryStatDir = "/".join(statFilePrefix.split("/")[:-1])

statH = {}
for fileName in os.listdir(summaryStatDir):
    if fileName.startswith(statFilePrefix.split("/")[-1]):
        winIndex = int(fileName.split("_")[1].split(".")[0])
        statFile = open(summaryStatDir + "/" + fileName)
        header = True
        assert not statH.has_key(winIndex)
        statH[winIndex] = {}
        for line in statFile.xreadlines():
            line = line.strip().split("\t")
            if header:
                header = False
                statNames = line
            else:
                for i in range(len(line)):
                    stat = float(line[i])
                    if not statH[winIndex].has_key(statNames[i]):
                        statH[winIndex][statNames[i]] = []
                    statH[winIndex][statNames[i]].append(stat)
        statFile.close()

winIndices = sorted(statH)
lenH = {}
for winIndex in winIndices:
    for statName in statH[winIndex].keys():
        lenH[len(statH[winIndex][statName])] = winIndex
if len(lenH) != 1:
    sys.exit("Not all windows have the same number of simulations (%s). AAAAAAAAAAAAARRRRRRRRRRRGGGGGGHHHHHHHHHHH!!!" %(lenH))

header = []
for j in range(len(statNames)):
    for k in range(len(statH)):
        header.append(statNames[j] + "_win%s" %(k))
print "\t".join(header)

for i in range(lenH.keys()[0]):#ith simulation
    outline = []
    for j in range(len(statNames)):#jth statistic
        statVec = []
        for k in range(len(statH)):#kth window
            statVec.append(statH[k][statNames[j]][i])
        minVal = min(statVec)
        if minVal < 0:
            statVec = [x-minVal for x in statVec]
        normStatVec = []
        #statMean = sum(statVec)/float(len(statVec))
        statSum = float(sum(statVec))
        if statSum == 0:
            normStatVec = [1.0/len(statVec)]*len(statVec)
        else:
            for k in range(len(statVec)):
                normStatVec.append(statVec[k]/statSum)
        outline.extend(normStatVec)
    print "\t".join([str(x) for x in outline])
