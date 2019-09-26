"""
================================================
My own Gaussion Mixture Model for VQSR.
================================================

Author: Shujia Huang
Date: 2017-08-02

"""
import re

import numpy as np
from sklearn.metrics import roc_curve

from basevar.log import logger
from basevar.io.openfile import Open

from basevar.caller.vqsr import variant_recalibrator_argument_collection as VRAC
from basevar.caller.vqsr import variant_datum as vd
from basevar.caller.vqsr import vcfutils


class VariantDataManager(object):

    def __init__(self, data=None):
        self.VRAC = VRAC.VariantRecalibratorArgumentCollection()
        self.annotation_mean = None
        self.annotation_STD = None
        # self.anno_texts = [['QUAL', 'Float', 'Raw variant quality before VQSR process'],
        #                   ['DP', 'Integer', 'Total depth of this variant'],
        #                   ['FS', 'Float', 'Phred-scaled p-value using '
        #                                   'Fisher\'s exact test to detect strand bias'],
        #                   ['Indel_SP', 'Integer', 'Indel species around this position.'
        #                                           'The less the better.'],
        #                   ['Indel_TOT', 'Integer', 'Number of Indel around this position.'
        #                                            'The less the better.']]

        self.data = []  # list <VariantDatum>
        if data:  # data is not None
            if not isinstance(data[0], vd.VariantDatum):
                raise ValueError('[ERROR] The data type should be "VariantDatum" in VariantDataManager(), '
                                 'but found %s' % type(data[0]))
            self.data = data
            for i, d in enumerate(self.data):
                self.data[i].annotations = np.array(self.data[i].annotations)

    def set_data(self, data):

        if not isinstance(data[0], vd.VariantDatum):
            raise ValueError('[ERROR] The data type should be "VariantDatum" in VariantDataManager(),'
                             'but found %s' % type(data[0]))
        self.data = data
        for i, d in enumerate(self.data):
            self.data[i].annotations = np.array(d.annotations)

    def normalization(self):
        # data normalization

        data = np.array([d.annotations for d in self.data], dtype=float)
        mean = data.mean(axis=0)
        self.annotation_mean = mean

        std = data.std(axis=0)
        self.annotation_STD = std

        # foundZeroVarianceAnnotation
        if any(std < 1e-5):
            raise ValueError('[ERROR] Found annotations with zero variance. '
                             'They must be excluded before proceeding.')

        # Each data now is (x - mean)/std
        for i, d in enumerate(data):

            self.data[i].annotations = (d - mean) / std

            # trim data by standard deviation threshold and mark failing data 
            # for exclusion later
            self.data[i].failing_STD_threshold = False
            if any(np.abs(self.data[i].annotations) > self.VRAC.STD_THRESHOLD):
                self.data[i].failing_STD_threshold = True

    def get_training_data(self):

        training_data = [d for d in self.data if ((not d.failing_STD_threshold) and d.at_training_site)]
        logger.info(('Training with %d variants after standard '
                     'deviation thresholding.\n' % len(training_data)))

        if len(training_data) < self.VRAC.MIN_NUM_BAD_VARIANTS:
            logger.warning('Training with very few variant sites! '
                           'Please check the model report and make '
                           'sure the quality of the model is reliable.')

        if len(training_data) > self.VRAC.MAX_NUM_TRAINING_DATA:
            logger.warning('Very large training set detected. '
                           'Downsampling to %d training variants.\n' %
                           self.VRAC.MAX_NUM_TRAINING_DATA)

            np.random.shuffle(training_data)  # Random shuffling
            return list([training_data[i] for i in range(self.VRAC.MAX_NUM_TRAINING_DATA)])

        return training_data

    def select_worst_variants(self, bad_lod):

        training_data = []
        for i, d in enumerate(self.data):
            if (d.lod < bad_lod) and (not d.failing_STD_threshold):
                training_data.append(d)
                # I do need: i order to be the same as self.data
                self.data[i].at_anti_training_site = True

        logger.info('Training with worst %d scoring variants --> variants with LOD < %.2f.\n' %
                    (len(training_data), bad_lod))

        if len(training_data) > self.VRAC.MAX_NUM_TRAINING_DATA:
            logger.warning('Very large training set detected.'
                           'Downsampling to %d training variants.\n' %
                           self.VRAC.MAX_NUM_TRAINING_DATA)

            np.random.shuffle(training_data)  # Random shuffling
            return list(training_data[i] for i in range(self.VRAC.MAX_NUM_TRAINING_DATA))

        return training_data

    def calculate_worst_lod_cutoff(self):

        lod_threshold, lod_cum = None, []
        if len(self.data) > 0:
            lod_dist = np.array([[d.at_training_site, d.lod]
                                 for d in self.data if (not d.failing_STD_threshold)])

            # I just use the 'roc_curve' function to calculate the worst 
            # LOD threshold, not use it to draw ROC curve And 'roc_curve' 
            # function will output the increse order, so that I don't 
            # have to sort it again

            _, tpr, thresholds = roc_curve(lod_dist[:, 0], lod_dist[:, 1])
            lod_cum = [[thresholds[i], 1.0 - r] for i, r in enumerate(tpr)]

            for i, r in enumerate(tpr):
                if r > 1.0 - self.VRAC.POSITIVE_TO_NEGATIVE_RATE:
                    lod_threshold = round(thresholds[i])
                    break

        return lod_threshold, np.array(lod_cum)


def load_training_site_from_VCF(vcf_file):
    """
    Just record the training site positions
    """
    logger.info('Loading Training site from VCF.\n')
    cdef int n = 0

    data_set = set()
    with Open(vcf_file, 'r') as I:
        for line in I:
            n += 1
            if n % 100000 == 0:
                logger.info("Loading lines %d" % n)

            if re.search(r'^#', line):
                continue

            col = line.strip().split()
            data_set.add(col[0] + ':' + col[1])  # just get the positions

    logger.info('[INFO] Finish loading training set %d lines.' % n)

    return data_set

def load_data_set(vcf_infile, training_set):
    if len(training_set) == 0:
        raise ValueError('[ERROR] No Training Data found')

    I = Open(vcf_infile, 'r')
    logger.info('Loading data set from VCF %s' % vcf_infile)

    n, data, h_info = 0, [], vcfutils.Header()
    with Open(vcf_infile, 'r') as I:
        for line in I:
            # VCF format
            n += 1
            if n % 100 == 0:
                logger.info('Loading lines %d' % n)

            # Record the header information
            if re.search(r'^#', line):
                h_info.record(line.strip())
                continue

            col = line.strip().split()
            if col[3] in ['N', 'n']:
                continue

            # qual = float(col[5])

            qd = re.search(r';?QD=([^;]+)', col[7])
            fs = re.search(r';?FS=([^;]+)', col[7])
            base_q_ranksum = re.search(r';?BaseQRankSum=([^;]+)', col[7])
            sor = re.search(r';?SOR=([^;]+)', col[7])
            mq_ranksum = re.search(r';?MQRankSum=([^;]+)', col[7])
            read_pos_ranksum = re.search(r';?ReadPosRankSum=([^;]+)', col[7])

            if any([not qd, not fs, not base_q_ranksum, not sor, not mq_ranksum, not read_pos_ranksum]):
                continue

            qd = round(float(qd.group(1)), 3) if qd.group(1) != "nan" else 10000
            fs = round(float(fs.group(1)), 3) if fs.group(1) != "nan" else 10000
            base_q_ranksum = round(float(base_q_ranksum.group(1)), 3) if base_q_ranksum.group(1) != "nan" else 10000
            sor = round(float(sor.group(1)), 3) if sor.group(1) != "nan" else 10000
            mq_ranksum = round(float(mq_ranksum.group(1)), 3) if mq_ranksum.group(1) != "nan" else 10000
            read_pos_ranksum = round(float(read_pos_ranksum.group(1)), 3) if read_pos_ranksum.group(
                1) != "nan" else 10000

            if fs >= 10000.0:
                fs = 10000.0

            datum = vd.VariantDatum()
            datum.raw_annotations = dict(QD=qd,
                                         FS=fs,
                                         BaseQRanksum=base_q_ranksum,
                                         SOR=sor,
                                         MQRanksum=mq_ranksum,
                                         ReadPosRankSum=read_pos_ranksum)
            datum.annotations = [qd, fs, base_q_ranksum, sor, mq_ranksum, read_pos_ranksum]

            datum.variant_order = col[0] + ':' + col[1]
            if datum.variant_order in training_set:
                datum.at_training_site = True

            data.append(datum)

    logger.info('Finish loading data set %d lines.' % n)
    return h_info, np.array(data)
