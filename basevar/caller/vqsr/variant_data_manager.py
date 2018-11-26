"""
================================================
My own Gaussion Mixture Model for SV genotyping.
Learn form scikit-learn
================================================

Author: Shujia Huang
Date: 2017-08-02

"""
import sys
import re
import time

import numpy as np
from sklearn.metrics import roc_curve

from . import variant_recalibrator_argument_collection as VRAC
from . import variant_datum as vd
from .. import vcfutils
from .. import utils


class VariantDataManager(object):

    def __init__(self, data=None):
        self.VRAC = VRAC.VariantRecalibratorArgumentCollection()
        self.annotationMean = None
        self.annotationSTD  = None
        self.annoTexts      = [['QUAL', 'Float', 'Raw variant quality before VQSR process'],
                               ['DP', 'Integer', 'Total depth of this variant'],
                               ['FS', 'Float', 'Phred-scaled p-value using '
                                               'Fisher\'s exact test to detect strand bias'],
                               ['Indel_SP', 'Integer', 'Indel species around this position.'
                                                       'The less the better.'],
                               ['Indel_TOT', 'Integer', 'Number of Indel around this position.'
                                                       'The less the better.']]

        self.data = [] # list <VariantDatum>
        if data: # data is not None
            if not isinstance(data[0],vd.VariantDatum):
                raise ValueError('[ERROR] The data type should be '
                                 '"VariantDatum" in VariantDataMa-'
                                 'nager(),but found %s'% str(type(data[0])))
            self.data = data
            for i, d in enumerate(self.data):
                self.data[i].annotations = np.array(self.data[i].annotations)

    def SetData(self, data):

        if not isinstance(data[0], vd.VariantDatum):
            raise ValueError('[ERROR] The data type should be "VariantDatum" '
                             'in VariantDataManager(),but found %s' %
                             str(type(data[0])))
        self.data = data
        for i, d in enumerate(self.data):
            self.data[i].annotations = np.array(d.annotations)

    def NormalizeData(self):

        data = np.array([d.annotations for d in self.data], dtype=float)
        mean = data.mean(axis=0)
        self.annotationMean = mean

        std  = data.std(axis=0)
        self.annotationSTD  = std

        # foundZeroVarianceAnnotation
        if any(std < 1e-5):
            raise ValueError('[ERROR] Found annotations with zero variance. '
                             'They must be excluded before proceeding.')

        # Each data now is (x - mean)/std
        for i, d in enumerate(data):

            self.data[i].annotations = (d - mean) / std
            # trim data by standard deviation threshold and mark failing data 
            # for exclusion later
            self.data[i].failingSTDThreshold = False
            if any(np.abs(self.data[i].annotations) > self.VRAC.STD_THRESHOLD):
                self.data[i].failingSTDThreshold = True

    def GetTrainingData(self):

        trainingData = [d for d in self.data if ((not d.failingSTDThreshold)
                                                 and d.atTrainingSite)]
        sys.stderr.write(('[INFO] Training with %d variants after standard '
                          'deviation thresholding.\n' % len(trainingData)))

        if len(trainingData) < self.VRAC.MIN_NUM_BAD_VARIANTS:
            sys.stderr.write(('[WARNING] Training with very few variant '
                              'sites! Please check the model reporting '
                              'PDF to ensure the quality of the model is '
                              'reliable.\n'))

        if len(trainingData) > self.VRAC.MAX_NUM_TRAINING_DATA:
            sys.stderr.write(('[WARING] Very large training set detected. '
                              'Downsampling to %d training variants.\n' %
                              self.VRAC.MAX_NUM_TRAINING_DATA))

            np.random.shuffle(trainingData) # Random shuffling
            return list(trainingData[i] 
                        for i in range(self.VRAC.MAX_NUM_TRAINING_DATA))

        return trainingData 

    def SelectWorstVariants(self, badLod):

        trainingData = []
        for i,d in enumerate(self.data):
            if(d.lod < badLod) and (not d.failingSTDThreshold):
                trainingData.append(d)
                # I do need: i order to be the same as self.data
                self.data[i].atAntiTrainingSite = True

        sys.stderr.write('[INFO] Training with worst %d scoring variants '
                         '--> variants with LOD < %.2f.\n' %
                         (len(trainingData), badLod))

        if len(trainingData) > self.VRAC.MAX_NUM_TRAINING_DATA:
            sys.stderr.write('[WARING] Very large training set detected.'
                             'Downsampling to %d training variants.\n' %
                             self.VRAC.MAX_NUM_TRAINING_DATA)

            np.random.shuffle(trainingData) # Random shuffling
            return list(trainingData[i] for i in range(self.VRAC.MAX_NUM_TRAINING_DATA))

        return trainingData

    def CalculateWorstLodCutoff(self):

        lodThreshold, lodCum = None, []
        if len(self.data) > 0:
            lodDist = np.array([[d.atTrainingSite, d.lod]
                                for d in self.data if(not d.failingSTDThreshold)])

            # I just use the 'roc_curve' function to calculate the worst 
            # LOD threshold, not use it to draw ROC curve And 'roc_curve' 
            # function will output the increse order, so that I don't 
            # have to sort it again
            _, tpr, thresholds = roc_curve(lodDist[:,0], lodDist[:,1]) 
            lodCum = [[thresholds[i], 1.0 - r] for i, r in enumerate(tpr)]

            for i, r in enumerate(tpr):
                if r > 1.0 - self.VRAC.POSITIVE_TO_NEGATIVE_RATE: 
                    lodThreshold = round(thresholds[i])
                    break

        return lodThreshold, np.array(lodCum)

def LoadTrainingSiteFromVCF(vcffile):
    """
    Just record the training site positions
    """
    I = utils.Open(vcffile, 'r')
    sys.stderr.write('\n[INFO] Loading Training site from VCF %s\n' % time.asctime())
    n, dataSet =0, set()
    for line in I:
        n += 1
        if n % 100000 == 0:
            sys.stderr.write('** Loading lines %d %s\n' % (n, time.asctime()))

        if re.search(r'^#', line):
            continue

        col = line.strip().split()
        dataSet.add(col[0] + ':' + col[1])  # just get the positions

    I.close()
    sys.stderr.write('[INFO] Finish loading training set %d lines. %s\n' %
                     (n, time.asctime()))

    return dataSet

def LoadDataSet(vcfInfile, traningSet):
    """
    """

    if len(traningSet) == 0:
        raise ValueError('[ERROR] No Training Data found')

    I = utils.Open(vcfInfile, 'r')
    sys.stderr.write('\n[INFO] Loading data set from VCF %s\n' % time.asctime())

    n, data, hInfo = 0, [], vcfutils.Header()
    for line in I: # VCF format
        n += 1
        if n % 100 == 0:
            sys.stderr.write('** Loading lines %d %s\n' % (n, time.asctime()))

        # Record the header information
        if re.search(r'^#', line):
            hInfo.record(line.strip())
            continue

        col = line.strip().split()
        if col[3] in ['N', 'n']:
            continue

        qual = float(col[5])

        dp = re.search(r';?CM_DP=([^;]+)', col[7])
        fs = re.search(r';?FS=([^;]+)', col[7])
        indel_sp = re.search(r';?Indel_SP=([^;]+)', col[7])
        indel_tot = re.search(r';?Indel_TOT=([^;]+)', col[7])
        if any([not dp, not fs, not indel_sp, not indel_tot]):
            continue

        dp = round(float(dp.group(1)), 2)
        fs = round(float(fs.group(1)), 3)
        indel_sp = float(indel_sp.group(1))
        indel_tot = float(indel_tot.group(1))

        if fs >= 10000.0:
            fs = 10000.0

        datum = vd.VariantDatum()
        datum.raw_annotations = dict(QUAL=qual,
                                     DP=dp,
                                     FS=fs,
                                     Indel_SP=indel_sp,
                                     Indel_TOT=indel_tot)

        datum.annotations = [qual, dp, fs, indel_sp, indel_tot]

        datum.variantOrder = col[0] + ':' + col[1]
        if datum.variantOrder in traningSet:
            datum.atTrainingSite = True

        data.append(datum)

    I.close()
    sys.stderr.write('[INFO] Finish loading data set %d lines. %s\n' %
                     (n, time.asctime()))

    return hInfo, np.array(data)

