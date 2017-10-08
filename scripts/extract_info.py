"""
================================================
Extract positive and negative training set data
================================================
"""
import os
import re
import sys

import time


def LoadDataSet(vcfInfile, traningSet):
    """
    """

    I = os.popen('gzip -dc %s' % vcfInfile) if vcfInfile.endswith('.gz') else open(vcfInfile)
    sys.stderr.write('\n[INFO] Loading data set from VCF %s\n' % time.asctime())

    n, data = 0, []
    for line in I: # VCF format
        n += 1
        if n % 100 == 0:
            sys.stderr.write('** Loading lines %d %s\n' % (n, time.asctime()))

        # Record the header information
        if re.search(r'^#', line):
            continue

        col = line.strip().split()
        qual = float(col[5])

        positive_set = re.search(r';?NEGATIVE_TRAIN_SITE=([^;]+)', col[7])
        negative_set = re.search(r';?POSITIVE_TRAIN_SITE=([^;]+)', col[7])

        if not positive_set or not negative_set:
            continue

        dp = re.search(r';?CM_DP=([^;]+)', col[7])
        qu = re.search(r';?QUAL=([^;]+)', col[7])
        fs = re.search(r';?FS=([^;]+)', col[7])

        if not dp or not fs:
            continue

        dp = round(float(dp.group(1)), 2)
        fs = round(float(fs.group(1)), 3)

        if fs >= 10000.0:
            fs = 10000.0

        datum.variantOrder = col[0] + ':' + col[1]
        if datum.variantOrder in traningSet:
            datum.atTrainingSite = True

        data.append(datum)

    I.close()

    sys.stderr.write('[INFO] Finish loading data set %d lines. %s\n' %
                     (n, time.asctime()))

    return hInfo, np.array(data)