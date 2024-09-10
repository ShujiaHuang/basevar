"""
=========================================
Variant quality score recalibrator (VQSR)
=========================================
Author: Shujia Huang & Siyang Liu
Date  : 2014-05-23 11:21:53
"""
import re
import sys
import time

from . import variant_data_manager as vdm
from . import variant_recalibrator as vror

from .. import utils


def main(opt):
    # Just record the sites of training data
    traningSet = vdm.LoadTrainingSiteFromVCF(opt.trainData)

    # Identify the traning sites
    hInfo, dataSet = vdm.LoadDataSet(opt.vcfInfile, traningSet)

    # init VariantRecalibrator object
    vr = vror.VariantRecalibrator()

    # Traning modul and calculate the VQ for all dataSet
    vr.OnTraversalDone(dataSet)
    # vr.VisualizationLodVStrainingSet(opt.figure + '.BadLodSelectInTraining')

    # For Record the Annnotations' values
    for d in vr.dataManager.annoTexts:
        hInfo.add('INFO', d[0], 1, d[1], d[2])

    # Outputting the result as VCF format
    hInfo.add('INFO', 'VQ', 1, 'Float', 'Variant Quality')
    hInfo.add('INFO', 'CU', 1, 'String', 'The annotation which was the worst '
                                         'performing in the Gaussian mixture modul, likely the reason why '
                                         'the variant was filtered out. It\'s the same tag as <culprit> '
                                         'in GATK')
    hInfo.add('INFO', 'NEGATIVE_TRAIN_SITE', 0, 'Flag',
              'This variant was used to build the negative training set of bad variants')
    hInfo.add('INFO', 'POSITIVE_TRAIN_SITE', 0, 'Flag',
              'This variant was used to build the positive training set of good variants')

    culprit, good, tot = {}, {}, 0.0
    annoTexts = [d[0] for d in vr.dataManager.annoTexts]

    for k, h in sorted(hInfo.header.items(), key=lambda d: d[0]):
        print (h)

    sys.stderr.write('\n[INFO] Outputting %s ...\n' % time.asctime())
    I = utils.Open(opt.vcfInfile, 'r')

    n, j, monitor = 0, 0, True
    for line in I:
        n += 1
        if n % 100000 == 0:
            sys.stderr.write('** Output lines %d %s\n' % (n, time.asctime()))

        if line.startswith('#'): continue

        col = line.strip().split()
        if col[3] in ['N', 'n']:
            continue

        dp = re.search(r';?CM_DP=([^;]+)', col[7])
        fs = re.search(r';?FS=([^;]+)', col[7])
        indel_sp = re.search(r';?Indel_SP=([^;]+)', col[7])
        indel_tot = re.search(r';?Indel_TOT=([^;]+)', col[7])
        if any([not dp, not fs, not indel_sp, not indel_tot]):
            continue

        order = col[0] + ':' + col[1]
        d = dataSet[j]
        j += 1  # Increase the index of dataSet for the next cycle
        if d.variantOrder != order:
            raise ValueError('[BUG] The order(%s) must be the same as '
                             'dataSet(%s)' % (order, d.variantOrder))

        # Deal with the INFO line
        vcfinfo = {}
        for info in col[7].split(';'):
            k = info.split('=')[0]

            if monitor and k in vcfinfo:
                monitor = False
                sys.stderr.write(('[WARNING] The tag: %s double hits in '
                                  'the INFO column at %s\n' %
                                  (k, opt.vcfInfile)))
            vcfinfo[k] = info

        tot += 1.0  # Record For summary
        culprit[annoTexts[d.worstAnnotation]] = culprit.get(
            annoTexts[d.worstAnnotation], 0.0) + 1.0  # For summary

        d.lod = round(d.lod * 10, 2)
        for lod in [0, 1, 2, 3, 4, 5, 10, 20, 25, 30, 35, 40, 45, 50]:
            if d.lod >= lod:
                good[lod] = good.get(lod, 0.0) + 1.0

        if d.atTrainingSite:
            vcfinfo['POSITIVE_TRAIN_SITE'] = 'POSITIVE_TRAIN_SITE'

        if d.atAntiTrainingSite:
            vcfinfo['NEGATIVE_TRAIN_SITE'] = 'NEGATIVE_TRAIN_SITE'

        vcfinfo['VQ'] = 'VQ=' + str(d.lod)
        vcfinfo['CU'] = 'CU=' + annoTexts[d.worstAnnotation]
        for k, v in d.raw_annotations.items():
            if k not in vcfinfo:
                vcfinfo[k] = k + '=' + str('%.2f' % v)

        col[7] = ';'.join(sorted(vcfinfo.values()))
        if d.lod < 0:
            d.lod = 0  # QUAL: donot allow less than 0

        col[5] = str(d.lod)  # QUAL field should use phred scala
        print ('\t'.join(col))

    I.close()

    sys.stderr.write('[INFO] Finish Outputting %d lines. %s\n' % (n, time.asctime()))

    ## Output Summary
    sys.stderr.write('\n[Summmary] Here is the summary information:\n')
    for k, v in sorted(good.items(), key=lambda k: k[0]):
        sys.stderr.write(('  ** Variant Site score >= %d: %d\t%0.2f\n' %
                          (k, v, v * 100 / tot)))

    for k, v in sorted(culprit.items(), key=lambda k: k[0]):
        sys.stderr.write(('  ** Culprit by %s: %d\t%.2f\n' % (k, v, v * 100.0 / tot)))
