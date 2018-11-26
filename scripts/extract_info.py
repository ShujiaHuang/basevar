"""
================================================
Extract positive and negative training set data
================================================
"""
import re
import sys
import time

import pysam


def LoadDataSet(vcfInfile):
    """
    """
    sys.stderr.write('\n[INFO] Loading data set from VCF %s\n' % time.asctime())

    positive_set = re.search(r';?NEGATIVE_TRAIN_SITE=([^;]+)', col[7])

    dp = round(float(dp.group(1)), 2)

    datum.variantOrder = col[0] + ':' + col[1]
    if datum.variantOrder in traningSet:
        datum.atTrainingSite = True

    data.append(datum)

    I.close()

    sys.stderr.write('[INFO] Finish loading data set %d lines. %s\n' %
                     (n, time.asctime()))

    return hInfo, np.array(data)