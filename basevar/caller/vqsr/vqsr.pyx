"""
=========================================
Variant quality score recalibrator (VQSR)
=========================================
Author: Shujia Huang & Siyang Liu
Date  : 2014-05-23 11:21:53
"""
import re

from basevar.log import logger
from basevar.caller.vqsr import variant_data_manager as vdm
from basevar.caller.vqsr import variant_recalibrator as vror

from basevar.io.openfile import Open

def main(opt):
    # Just record the sites of training data
    training_set = vdm.load_training_site_from_VCF(opt.train_data)

    # Identify the training sites
    h_info, data_set = vdm.load_data_set(opt.vcf_infile, training_set)

    # init VariantRecalibrator object
    vr = vror.VariantRecalibrator()

    # Training model and calculate the VQ for all data_set
    vr.on_traversal_done(data_set)
    # vr.visualization_lod_VS_training_set(opt.figure + '.BadLodSelectInTraining')

    # For Record the Annnotations' values
    # for d in vr.data_manager.anno_texts:
    #     h_info.add('INFO', d[0], 1, d[1], d[2])

    # Outputting the result as VCF format
    h_info.add('INFO', 'VQSLOD', 1, 'Float', 'Variant quality calculate by VQSR')
    h_info.add('INFO', 'CU', 1, 'String',
               'The annotation which was the worst '
               'performing in the Gaussian mixture modul, likely the reason why '
               'the variant was filtered out. It\'s the same tag as <culprit> '
               'in GATK')
    h_info.add('INFO', 'NEGATIVE_TRAIN_SITE', 0, 'Flag',
               'This variant was used to build the negative training set of bad variants')
    h_info.add('INFO', 'POSITIVE_TRAIN_SITE', 0, 'Flag',
               'This variant was used to build the positive training set of good variants')

    # anno_texts = [d[0] for d in vr.data_manager.anno_texts]
    culprit, good, tot = {}, {}, 0.0
    anno_texts = ['QD', 'FS', 'BaseQRankSum', 'SOR', 'MQRankSum', 'ReadPosRankSum']

    for k, h in sorted(h_info.header.items(), key=lambda d: d[0]):
        print (h)

    logger.info("[INFO] Outputting %s ...")
    I = Open(opt.vcf_infile, 'r')

    n, j, monitor = 0, 0, True
    for line in I:
        n += 1
        if n % 100000 == 0:
            logger.info("** Output lines %d." % n)

        if line.startswith('#'):
            continue

        col = line.strip().split()
        if col[3] in ['N', 'n']:
            continue

        qd = re.search(r';?QD=([^;]+)', col[7])
        fs = re.search(r';?FS=([^;]+)', col[7])
        base_q_ranksum = re.search(r';?BaseQRankSum=([^;]+)', col[7])
        sor = re.search(r';?SOR=([^;]+)', col[7])
        mq_ranksum = re.search(r';?MQRankSum=([^;]+)', col[7])
        read_pos_ranksum = re.search(r';?ReadPosRankSum=([^;]+)', col[7])

        if any([not qd, not fs, not base_q_ranksum, not sor, not mq_ranksum, not read_pos_ranksum]):
            continue

        order = col[0] + ':' + col[1]
        d = data_set[j]
        j += 1  # increase the index of data_set for the next cycle.
        if d.variant_order != order:
            raise ValueError('[BUG] The order(%s) must be the same as '
                             'dataSet(%s)' % (order, d.variant_order))

        # Deal with the INFO line
        vcf_info = {}
        for info in col[7].split(';'):
            k = info.split('=')[0]

            if monitor and k in vcf_info:
                monitor = False
                logger.warning('The tag: %s double hits in the INFO column at %s.' %
                                (k, opt.vcf_infile))
            vcf_info[k] = info

        tot += 1.0  # Record For summary
        culprit[anno_texts[d.worst_annotation]] = culprit.get(
            anno_texts[d.worst_annotation], 0.0) + 1.0  # For summary

        d.lod = round(d.lod * 10, 2)
        for lod in [0, 1, 2, 3, 4, 5, 10, 20, 25, 30, 35, 40, 45, 50]:
            if d.lod >= lod:
                good[lod] = good.get(lod, 0.0) + 1.0

        if d.at_training_site:
            vcf_info['POSITIVE_TRAIN_SITE'] = 'POSITIVE_TRAIN_SITE'

        if d.at_anti_training_site:
            vcf_info['NEGATIVE_TRAIN_SITE'] = 'NEGATIVE_TRAIN_SITE'

        vcf_info['VQSLOD'] = 'VQSLOD=' + str(d.lod)
        vcf_info['CU'] = 'CU=' + anno_texts[d.worst_annotation]
        for k, v in d.raw_annotations.items():
            if k not in vcf_info:
                vcf_info[k] = k + '=' + str('%.2f' % v)

        col[7] = ';'.join(sorted(vcf_info.values()))
        if d.lod < 0:
            d.lod = 0  # QUAL: donot allow less than 0

        # col[5] = str(d.lod)
        print ('\t'.join(col))

    I.close()

    logger.info('Finish Outputting %d lines.\n' % n)

    ## Output Summary
    logger.info('[Summmary] Here is the summary information:')
    for k, v in sorted(good.items(), key=lambda k: k[0]):
        logger.info(('  ** Variant Site score >= %d: %d\t%0.2f' %
                     (k, v, v * 100 / tot)))

    for k, v in sorted(culprit.items(), key=lambda k: k[0]):
        logger.info(('  ** Culprit by %s: %d\t%.2f' % (k, v, v * 100.0 / tot)))
