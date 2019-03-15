"""
===========================================
===========================================
Author: Shujia Huang & Siyang Liu
Date  : 2014-05-20 08:50:06
"""
import sys
import numpy as np
from sklearn.mixture import GaussianMixture as GMM
from sklearn.utils.extmath import logsumexp

# My own class
from . import variant_datum as vd
from . import variant_recalibrator_argument_collection as VRAC


class VariantRecalibratorEngine(object):

    def __init__(self, vrac=None):

        self.VRAC = VRAC.VariantRecalibratorArgumentCollection()
        if vrac:
            self.VRAC = vrac

        self.MIN_PROB_CONVERGENCE = 2e-3
        self.MIN_ACCEPTABLE_LOD_SCORE = -2000.0

    def ClassifyData(self, dataSize):
        """
        Classify the data into TrainingSet, Cross-ValidationSet and TestSet. 
        Reture the data indexes

        Call in GenerateModel
        """
        trainSetSize = int(np.round(self.VRAC.TRAIN_SIZE_RATE * dataSize))
        cvSetSize = int(np.round(self.VRAC.CV_SIZE_RATE * dataSize))
        testSetSize = int(np.round(self.VRAC.TEST_SIZE_RATE * dataSize))

        # The index array of training data
        trainSetIdx = range(trainSetSize)
        # The index array of cross-validation data 
        cvSetIdx = range(trainSetSize, cvSetSize + trainSetSize)
        # The index array of Test data
        testSetIdx = range(cvSetSize + testSetSize, dataSize)

        return trainSetIdx, cvSetIdx, testSetIdx

    def GenerateModel(self, data, maxGaussians):

        if len(data) == 0:
            raise ValueError('[ERROR] No data found. The size is %d\n' % len(data))

        if not isinstance(data[0], vd.VariantDatum):
            raise ValueError('[ERROR] The data type should be "VariantDatum" '
                             'in GenerateModel() of class VariantRecalibrato-'
                             'rEngine(), but found %s\n' % str(type(data[0])))

        if maxGaussians <= 0:
            raise ValueError('[ERROR] maxGaussians must be a positive integer '
                             'but found: %d\n' % maxGaussians)

        gmms = [GMM(n_components=n + 1,
                    covariance_type='full',
                    tol=self.MIN_PROB_CONVERGENCE,
                    n_iter=self.VRAC.NITER,
                    n_init=self.VRAC.NINIT) for n in range(maxGaussians)]

        # gmms = [BGMM(n_components=n + 1,
        #              covariance_type='full',
        #              tol=self.MIN_PROB_CONVERGENCE,
        #              max_iter=self.VRAC.NITER,
        #              n_init=self.VRAC.NINIT) for n in range(maxGaussians)]

        trainingData = np.array([d.annotations for d in data])

        # np.random.shuffle(trainingData) # Random shuffling
        # trainSetIdx, cvSetIdx, testSetIdx = self.ClassifyData(len(trainingData))

        # find a best components for GMM model
        minBIC, bics = np.inf, []
        for g in gmms:
            sys.stderr.write('[INFO] Trying %d gaussian in GMM process '
                             'training ...\n' % g.n_components)

            g.fit(trainingData)

            bic = g.bic(trainingData)
            bics.append(bic)

            if bic == float('inf') or (bic < minBIC and g.converged_):
                bestgmm, minBIC = g, bic

            sys.stderr.write('  -- Converge infomation of training '
                             'process: %s\n' % g.converged_)

        sys.stderr.write('[INFO] All the BIC: %s\n' % bics)
        sys.stderr.write('[INFO] Model Training Done. And take the model '
                         'with %d gaussiones which with BIC %f.\n' %
                         (len(bestgmm.means_), minBIC))

        return bestgmm

    def EvaluateData(self, data, gmm, evaluateContrastively=False):

        if not isinstance(data[0], vd.VariantDatum):
            raise ValueError('[ERROR] The data type should be "VariantDatum" '
                             'in EvaluateData() of class VariantRecalibrator-'
                             'Engine(), but found %s\n' % str(type(data[0])))

        sys.stderr.write('[INFO] Evaluating full set of %d variants ...\n' % len(data))

        for i, _ in enumerate(data):

            # log likelihood and the base is 10
            thisLod = gmm.score(data[i].annotations[np.newaxis, :]) / np.log(10)
            if np.math.isnan(thisLod):
                gmm.converged_ = False
                return

            if evaluateContrastively:
                # data[i].lod must has been assigned by good model or something 
                # like that
                # contrastive evaluation: (prior + positive model - negative model)
                data[i].lod = data[i].prior + data[i].lod - thisLod
                if thisLod == float('inf'):
                    data[i].lod = self.MIN_ACCEPTABLE_LOD_SCORE * \
                                  (1.0 + np.random.rand(1)[0])
            else:
                # positive model only so set the lod and return 
                data[i].lod = thisLod

        return self

    def CalculateWorstPerformingAnnotation(self, data, goodModel, badModel):

        for i, d in enumerate(data):
            probDiff = [self.EvaluateDatumInOneDimension(goodModel, d, k) -
                        self.EvaluateDatumInOneDimension(badModel, d, k)
                        for k in range(len(d.annotations))]

            # Get the index of the worst annotations
            data[i].worstAnnotation = np.argsort(probDiff)[0]

        return self

    def EvaluateDatumInOneDimension(self, gmm, datum, iii):

        # pVarInGaussianLogE = [
        #         np.log(w) + NormalDistributionLoge(
        #                 gmm.means_[k][iii],
        #                 gmm.covars_[k][iii][iii],
        #                datum.annotations[iii])
        #        for k, w in enumerate(gmm.weights_)
        # ]
        pVarInGaussianLogE = [
            np.log(w) + NormalDistributionLoge(
                gmm.means_[k][iii],
                gmm.covars_[k][iii][iii],
                datum.annotations[iii])
            for k, w in enumerate(gmm.weights_)
        ]

        # np.log10(Sum(pi_k * p(v|n,k)))
        return logsumexp(np.array(pVarInGaussianLogE)) / np.log(10)


def NormalDistributionLoge(mu, sigma, x):
    if sigma <= 0:
        raise ValueError('[ERROR] sd: Standard deviation of normal must '
                         'be > 0 but found: %f\n' % sigma)
    if (mu == float('inf') or mu == float('-inf') or
            sigma == float('inf') or sigma == float('-inf') or
            x == float('inf') or x == float('-inf')):
        raise ValueError('[ERROR] mean, sd, or, x: Normal parameters must '
                         'be well formatted (non-INF, non-NAN)')

    a = -1.0 * (np.log(sigma) + 0.5 * np.log(2 * np.pi))
    b = -0.5 * ((x - mu) / sigma) ** 2

    return a + b  # The Natural log
