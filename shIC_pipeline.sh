#!/bin/bash
# This file describes the pipeline for training and running S/HIC, along with an executable example.
# All of the commands in the examples below are executed from the shIC_code/ directory.

# STEP ONE: Build C tools needed for calculating summary statistics
# To build niceStats and pgStatsBed, which we use to calculate population genetic summary statistics
# on simulated and real data, respectively, simply navigate to the S/HIC project folder and run the
# following command:
make

# STEP TWO: generate training feature vectors from simulations:
# The bash loop below requires a directory with ms-formatted simulation output files 
# (trainingSetsTennessenEuro in the example below). If we are going to use 11 subwindows for our 
# classifier, this directory must have 23 files: 11 files with simulated hard sweeps whose file 
# names end with "Hard_0.msout", "Hard_1.msout" and so on, all the way up to "Hard_10.msout"
# The number in the file corresponds to the location of the sweep in each simulation: for the file
# ending with "Hard_0.msout" the sweep occurs in the middle of the first window, for "Hard_1.msout"
# the sweep is in the middle of the second window, and so on. Similarly, there are 11 files with
# simulated soft sweeps, ending with "Soft_0.msout", etc. Finally, there is one file with neutral
# simulations, whose name ends with "Neut.msout"
# The bash loop below splits these simulations into 11 subwindows (using
# splitMsOutputIntoWindows.py), calculates summary statistics in each of these subwindows
# (via niceStats), and then combines the statistics from the 11 subwindows into one feature
# vector for each training example (via combineWinsIntoFeatureVec.py). The user can modify this 
# loop to use training data from a different directory <dir> by replacing each instance of 
# "trainingSetsTennessenEuro" in the code below with <dir>. Running this loop will then create a 
# directory, <dir>FeatureVecs, that contains feature vectors corresponding to each simulated 
# training instance. This step may take some time (on the order of 30 minutes for fairly large 
# simulated chromosomes with 1000 training examples in each file.

# loop for going from ms output to training feature vecs:	
mkdir trainingSetsTennessenEuroFeatureVecs
for file in `ls trainingSetsTennessenEuro`;
do
	gunzip trainingSetsTennessenEuro/$file*
	file=${file/\.msout.gz/.msout}
	#split the training examples into 11 equally-sized adjacent windows
	python splitMsOutputIntoWindows.py trainingSetsTennessenEuro/$file 11 tmpMSFileWins;
        
	#run niceStats on each window
        mkdir tmpStatDir
	for winFile in `ls tmpMSFileWins*`;
        do
            cat $winFile | ./niceStats > tmpStatDir/$winFile.stats
            rm $winFile
        done
	
	#combine stats from individual windows into one set of feature vectors
	python combineWinsIntoFeatureVec.py tmpStatDir/tmpMSFileWins > trainingSetsTennessenEuroFeatureVecs/${file/msout/fvec}
        rm -rf tmpStatDir/
	gzip trainingSetsTennessenEuro/$file*
done

# STEP THREE: combine feature vectors for hard and soft sweeps of varying locations, and neutral sweeps into 5 training sets
# this step is performed by formatTrainingSets.py
# Arguments:
#	1) the file containing feature vectors for our neutral training examples
#	2) the prefix of all files containing feature vectors for hard sweep examples
#	3) the prefix of all files containing feature vectors for soft sweep examples
#	4) the index of the central subwindow (5, if we are using 11 subwindows)
#	5) the output directory, which will have 5 files (one for each training set) after running
mkdir combinedTrainingSetsTennessenEuro
python formatTrainingSets.py trainingSetsTennessenEuroFeatureVecs/tennessenEuroNeut.fvec trainingSetsTennessenEuroFeatureVecs/tennessenEuroSoft trainingSetsTennessenEuroFeatureVecs/tennessenEuroHard 5 combinedTrainingSetsTennessenEuro/

# STEP FOUR: train the classifier
# Next, we have to train the classifier, which is saved as a "pickle" using scikit-learn's joblib library
# Beware that pickles can contain executable python code, so if your pickles are tampered with there is the potential
# for execution of malicious code. Make sure your pickles cannot be edited by anyone else!
# For more information on pickles and joblib, see http://scikit-learn.org/stable/modules/model_persistence.html
# The joblib pickles can consist of many files, so it is best to save them each in their own directory, as done below.
mkdir classifiers
mkdir classifiers/tennessenEuro

# trainClassifier.py will train the Extra-Trees classifier using scikit-learn, and save it as a pickle using joblib
# Arguments:
#	1) directory with the 5 training sets, which should be named hard.fvec, soft.fvec, linkedHard.fvec, linkedSoft.fvec, and neut.fvec
#	2) the path/name of the pickle where the classifier will be stored
#	Arguments 3+): the remaining arguments are the names of statistics you wish to be used by the classifier.
#	They must match the names in the headers of the training feature vector, but omit the "_win<i>" suffix.
#	Simply supply the word "all" (without quotes) as argument 3 if you wish to use all statistics present in
#	the feature vector files.
python trainClassifier.py combinedTrainingSetsTennessenEuro/ classifiers/tennessenEuro/tennessenEuroAutosomal.p all
# Note: this step can done using two pre-computed feature vector directories included in shIC_code:
# combinedTrainingSetsTennessenEuro and combinedTrainingSetsEquib

# STEP FIVE: calculate summary statistics in genomic windows
# Before we classify genomic data, we must calculate summary statistics and use them to create our feature vectors.
# First, we have to calculate summary statistics, which can be done with pgStatsBed.
# This step must be performed on only one chromosome (or a portion of a chromosome) at a time.

#first, we have to combine the two parts of the sample population's fasta file (which was too large for github), then unzip
cat test2LCTpop_CEU.fa.gz.partaa test2LCTpop_CEU.fa.gz.partab > test2LCTpop_CEU.fa.gz
rm test2LCTpop_CEU.fa.gz.partaa test2LCTpop_CEU.fa.gz.partab
gunzip test2LCTpop_CEU.fa.gz

# pgStatsBed calculates summary statistics in individual subwindows
# Arguments:
#	1) the path to the ingroup fasta file, which contains the fasta-formatted sequences of each individual
#		Note: sites that you wish to omit from calculations (e.g. repeats) should be Ns in each individual
#	2) the path to the ancestral fasta file, which contains the inferred ancestral state (used for Fay and Wu's H and theta-H)
#		Note: simply use a fasta file with the reference sequence if you are not using Fay and Wu's statistics
#	3) the path to a bed-formatted file containing the coordinates of subwindows in which to calculate summary stats
#		Note: the bed format specification can be found here: http://genome.ucsc.edu/FAQ/FAQformat.html#format1
#	4) the maximum fraction of sites that can be missing from each subwindow; subwindows not meeting this criterion are skipped
./pgStatsBed test2LCTpop_CEU.fa test2LCTanc.fa chr2_200kb_134M_140M.bed 0.25 > chr2_200kb_134M_140M_CEU.ss

# STEP SIX: combine our summary statistics from individual genomic subwindows into feature vectors.
# bedSummaryStatsToFeatureVec.py generates feature vectors from summary statistic output
# i.e. once statistics have been calculated in subwindows, this file walks along the chromosome
# examining runs of windows, creating a feature vector to represent the central window as described in
# the paper.
# Arguments:
#	1) the path to the file with summary statistic values we wish to convert to feature vectors
#	2) the minimum number of unmasked sites each subwindow must have in order to be included
#	3) the minimum value of nucleotide diversity (pi) across a window for it to be included (value does not matter if pi is absent from file)
#	4) the number of subwindows (usually 11)
#	5) the stepSize between windows to classify (set to the window size to classify adjacent windows)
python bedSummaryStatsToFeatureVec.py chr2_200kb_134M_140M_CEU.ss 50000 0 11 200000 > chr2_200kb_134M_140M_CEU.fvec

# Warning: if you are using your own summary statistics output into bedSummaryStatsToFeatureVec.py, 
# note that currently, this script may change the values of some of your statistics depending on the 
# names given to them in the header line. This is because bedSummaryStatsToFeatureVec.py  currently 
# alters the values of some summary statistics generated by pgStatsBed, listed below with their header 
# names:
#	pi: pgStatsBed currently calculates total nucleotide diversity across the entire window, we divide by the number of unmasked sites in the window
#	thetaH: we also divide Fay and Wu's thetaH from the number of unmasked sites in the window
#	segSites: rather than Watterson's theta, pgStatsBed counts the number of segregating sites, which we divide by the total number of unmasked sites
#		this will yield an identical feature vector to Watterson's theta provided the sample size is constant across sites
#	HapCount: the number of distinct haplotypes is also divided by the number of unmasked sites prior to generating feature vectors
#	fayWuH: if Fay and Wu's H is included, it is recalculated by subtracting pi from thetaH (after each have been divided by the number of unmasked sites.
# If you are using as input a summary statistic file whose headers do not match any of these header names above, no values will be altered.

# STEP SEVEN: classify genomic data, which we can do using either of the following two commands:

# classifyWinsInChromosome.py classifies genomic data and print result in bed format
# Arguments:
#	1) the path to the classifier (pickle)
#	2) the path to file with feature vectors to classify
#	3) size of windows to classify (in kilo bases)
#	4) output format
# For the example below, the output class for each classified window is shown in the "Name" field (4th column).
# The resulting .bed file can be uploaded to the UCSC Genome Browser (http://genome.ucsc.edu/cgi-bin/hgGateway) as a custom track for visualization.
python classifyWinsInChromosome.py classifiers/tennessenEuro/tennessenEuroAutosomal.p chr2_200kb_134M_140M_CEU.fvec 200 bed > chr2_200kb_134M_140M_CEU.classifications.bed
# Note: for this example, we are examining a portion of chromosome 2 beginning at 134,000,000
# but the coordinates in the .bed, .ss, and .fvec files must begin with 0. So 134,000,000 must
# be added to each coordinate in chr2_200kb_134M_140M_CEU.classifications.bed to see the true
# coordinates, which can be done as follows:
awk 'BEGIN {count=0;} {if (count == 0) print $0; else printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2+134000000, $3+134000000, $4, $5, $6, $7+134000000, $8+134000000, $9; count += 1;}' chr2_200kb_134M_140M_CEU.classifications.bed > tmp.bed ; mv tmp.bed chr2_200kb_134M_140M_CEU.classifications.bed


# classifyWinsInChromosome.py can also classify genomic data and write posterior probabilities for each of the five classes in bedGraph format.
# 	The bedGraph format specification is available here: http://genome.ucsc.edu/FAQ/FAQformat.html#format1.8
# Here we change the output format from "bed" to "bedGraph", and provide an additional argument: 
# the desired prefix of the output files. This command will generate five different output files, 
# one for each class. The resulting .bedGraph files can be uploaded to the UCSC Genome Browser
# (http://genome.ucsc.edu/cgi-bin/hgGateway) as a custom track for visualization.
python classifyWinsInChromosome.py classifiers/tennessenEuro/tennessenEuroAutosomal.p chr2_200kb_134M_140M_CEU.fvec 200 bedGraph chr2_200kb_134M_140M_CEU.classProbs
# NOTE: again, for this example we need to add 134,000,000 to each coordinate:
for file in chr2_200kb_134M_140M_CEU.classProbs.*.bedGraph ; do awk 'BEGIN {count=0;} {if (count == 0) print $0; else printf "%s\t%s\t%s\t%s\n", $1, $2+134000000, $3+134000000, $4; count += 1;}' $file > tmp.bedGraph ; mv tmp.bedGraph $file ; done
