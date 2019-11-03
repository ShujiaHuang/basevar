"""
Author: Shujia Huang
Date  : 2014-05-23 11:21:53
"""
from basevar.log import logger
from basevar.caller.vqsr import variant_data_manager as vdm
from basevar.caller.vqsr import variant_recalibrator_engine as vre
from basevar.caller.vqsr import variant_recalibrator_argument_collection as VRAC


class VariantRecalibrator(object):

    def __init__(self):
        self.VRAC = VRAC.VariantRecalibratorArgumentCollection()
        self.data_manager = vdm.VariantDataManager()
        self.engine = vre.VariantRecalibratorEngine(self.VRAC)
        self.bad_lod_cutoff = None
        self.lod_cum_in_train = []

    def on_traversal_done(self, data):
        self.data_manager.set_data(data)
        self.data_manager.normalization()

        # Generate the positive model using the training data and evaluate 
        # each variant
        positive_training_data = self.data_manager.get_training_data()
        logger.info('\nTraining the goodModel ...')

        good_model = self.engine.generate_model(positive_training_data, self.VRAC.MAX_GAUSSIANS)

        logger.info('The converged information of goodModel is: %s.' % good_model.converged_)
        logger.info('The means of gaussion of goodModel is:\n%s.' % good_model.means_)

        self.engine.evaluate_data(self.data_manager.data, good_model, False)
        self.bad_lod_cutoff, self.lod_cum_in_train = self.data_manager.calculate_worst_lod_cutoff()

        # Generate the negative model using the worst performing data and 
        # evaluate each variant contrastively
        logger.info('\nTraining the badModel ...')
        negative_training_data = self.data_manager.select_worst_variants(self.bad_lod_cutoff)
        bad_model = self.engine.generate_model(
            negative_training_data,
            min(self.VRAC.MAX_GAUSSIANS_FOR_NEGATIVE_MODEL, self.VRAC.MAX_GAUSSIANS)
        )

        logger.info('The converged information of badModel is: %s.' % bad_model.converged_)
        logger.info('The means of gaussion of badModel is:\n%s.' % bad_model.means_)
        self.engine.evaluate_data(self.data_manager.data, bad_model, True)

        if (not good_model.converged_) or (not bad_model.converged_):
            raise ValueError('[ERROR] NaN LOD value assigned. Clustering '
                             'with these variants and these annotations is '
                             'unsafe. Please consider raising the number of '
                             'variants used to train the negative model or '
                             'lowering the maximum number of gaussians allowed '
                             'for use in the model.')

        # Find the VQSLOD cutoff values which correspond to the various 
        # tranches of calls requested by the user
        self.engine.calculate_worst_performing_annotation(self.data_manager.data, good_model, bad_model)

    def visualization_lod_VS_training_set(self, fig_name):
        import matplotlib.pyplot as plt

        fig = plt.figure()
        plt.title('LOD VS Positive training set', fontsize=14)
        plt.plot(self.lod_cum_in_train[:, 0], self.lod_cum_in_train[:, 1], 'r-')
        plt.scatter(self.lod_cum_in_train[:, 0], self.lod_cum_in_train[:, 1], c='r',
                    marker='.', linewidth=0, alpha=0.5)

        plt.plot([self.bad_lod_cutoff, self.bad_lod_cutoff], [0, 1], 'g--')
        plt.ylim(0, 1.0)
        plt.xlim(-10, 10)
        plt.xlabel('Variant score threshold for the bad model', fontsize=16)
        plt.ylabel('Rate of Positive->Negative', fontsize=16)

        fig.savefig(fig_name)
