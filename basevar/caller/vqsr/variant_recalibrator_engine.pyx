"""
Author: Shujia Huang & Siyang Liu
Date  : 2014-05-20 08:50:06
"""
import numpy as np
from scipy.misc import logsumexp
from sklearn.mixture import GaussianMixture

# My own class
from basevar.log import logger
from basevar.caller.vqsr import variant_datum as vd
from basevar.caller.vqsr import variant_recalibrator_argument_collection as VRAC


class VariantRecalibratorEngine(object):

    def __init__(self, vrac=None):

        self.VRAC = VRAC.VariantRecalibratorArgumentCollection()
        if vrac:
            self.VRAC = vrac

        self.MIN_PROB_CONVERGENCE = 2e-3
        self.MIN_ACCEPTABLE_LOD_SCORE = -2000.0

    def classify_data(self, data_size):
        """
        Classify the data into TrainingSet, Cross-ValidationSet and TestSet. 
        Reture the data indexes

        Call in GenerateModel
        """
        train_set_size = int(np.round(self.VRAC.TRAIN_SIZE_RATE * data_size))
        cv_set_size = int(np.round(self.VRAC.CV_SIZE_RATE * data_size))
        test_set_size = int(np.round(self.VRAC.TEST_SIZE_RATE * data_size))

        # The index array of training data
        train_set_idx = range(train_set_size)

        # The index array of cross-validation data
        cv_set_idx = range(train_set_size, cv_set_size + train_set_size)

        # The index array of Test data
        test_set_idx = range(cv_set_size + test_set_size, data_size)

        return train_set_idx, cv_set_idx, test_set_idx

    def generate_model(self, data, max_gaussians):

        if len(data) == 0:
            raise ValueError('[ERROR] No data found. The size is %d\n' % len(data))

        if not isinstance(data[0], vd.VariantDatum):
            raise ValueError('[ERROR] The data type should be "VariantDatum" '
                             'in GenerateModel() of class VariantRecalibrato-'
                             'rEngine(), but found %s\n' % str(type(data[0])))

        if max_gaussians <= 0:
            raise ValueError('[ERROR] maxGaussians must be a positive integer '
                             'but found: %d\n' % max_gaussians)

        gmms = [GaussianMixture(n_components=n + 1,
                                covariance_type='full',
                                tol=self.MIN_PROB_CONVERGENCE,
                                max_iter=self.VRAC.NITER,
                                n_init=self.VRAC.NINIT) for n in range(max_gaussians)]

        training_data = np.array([d.annotations for d in data])

        # find a best components for GMM model
        min_bic, bics = np.inf, []
        for g in gmms:
            logger.info('Trying %d gaussian in GMM process training ...' % g.n_components)

            g.fit(training_data)

            bic = g.bic(training_data)
            bics.append(bic)

            if bic == float('inf') or (bic < min_bic and g.converged_):
                best_gmm, min_bic = g, bic

            logger.info('  -- Converge information of training process: %s' % g.converged_)

        logger.info('[INFO] All the BIC: %s' % bics)
        logger.info('[INFO] Model Training Done. And take the model '
                    'with %d gaussiones which with BIC %f.\n' %
                    (len(best_gmm.means_), min_bic))

        return best_gmm

    def evaluate_data(self, data, gmm, evaluate_contrastively=False):

        if not isinstance(data[0], vd.VariantDatum):
            raise ValueError('[ERROR] The data type should be "VariantDatum" '
                             'in EvaluateData() of class VariantRecalibrator-'
                             'Engine(), but found %s' % str(type(data[0])))

        logger.info('Evaluating full set of %d variants ...' % len(data))

        for i, _ in enumerate(data):

            # log likelihood and the base is 10
            this_lod = gmm.score(data[i].annotations[np.newaxis, :]) / np.log(10)
            if np.math.isnan(this_lod):
                gmm.converged_ = False
                return

            if evaluate_contrastively:
                # data[i].lod must has been assigned by good model.
                # contrastive evaluation: (prior + positive model - negative model)
                data[i].lod = data[i].prior + data[i].lod - this_lod
                if this_lod == float('inf'):
                    data[i].lod = self.MIN_ACCEPTABLE_LOD_SCORE * (1.0 + np.random.rand(1)[0])
            else:
                # positive model only so set the lod and return 
                data[i].lod = this_lod

        return self

    def calculate_worst_performing_annotation(self, data, good_model, bad_model):

        for i, d in enumerate(data):
            prob_diff = [self.evaluate_datum_in_one_dimension(good_model, d, k) -
                         self.evaluate_datum_in_one_dimension(bad_model, d, k)
                         for k in range(len(d.annotations))]

            # Get the index of the worst annotations
            data[i].worst_annotation = np.argsort(prob_diff)[0]

        return self

    def evaluate_datum_in_one_dimension(self, gmm, datum, iii):

        p_var_in_gaussian_loge = [
            np.log(w) + normal_distribution_Loge(
                gmm.means_[k][iii],
                gmm.covariances_[k][iii][iii],  # gmm.covars_[k][iii][iii],
                datum.annotations[iii])
            for k, w in enumerate(gmm.weights_)
        ]

        # np.log10(Sum(pi_k * p(v|n,k)))
        return logsumexp(np.array(p_var_in_gaussian_loge)) / np.log(10)


def normal_distribution_Loge(mu, sigma, x):
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
