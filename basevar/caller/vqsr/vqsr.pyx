"""
=========================================
Variant quality score recalibrator (VQSR)
=========================================
Author: Shujia Huang & Siyang Liu
Date  : 2014-05-23 11:21:53
"""
import sys
import re

from basevar.log import logger
from basevar.caller.vqsr import variant_data_manager as vdm
from basevar.caller.vqsr import variant_recalibrator as vror

from basevar.io.openfile import Open
from basevar.io.BGZF.tabix import tabix_index

def run_VQSR(opt):
    # Just record the sites of training data
    training_set = vdm.load_training_site_from_VCF(opt.train_data)

    # Identify the training sites
    h_info, data_set = vdm.load_data_set(opt.vcf_infile, training_set)

    # init VariantRecalibrator object
    vr = vror.VariantRecalibrator()

    # Training model and calculate the VQ for all data_set
    vr.on_traversal_done(data_set)
    # vr.visualization_lod_VS_training_set('VQSR.Training.BadLodSelectInTraining.png')

    # Outputting the result as VCF format
    h_info.add('INFO', 'VQSLOD', 1, 'Float', 'Variant quality calculate by VQSR')
    h_info.add('INFO', 'CU', 1, 'String',
               'The annotation which was the worst performing in the Gaussian mixture module,'
               'likely the reason why the variant was filtered out.')
    h_info.add('INFO', 'NEGATIVE_TRAIN_SITE', 0, 'Flag',
               'This variant was used to build the negative training set of bad variants')
    h_info.add('INFO', 'POSITIVE_TRAIN_SITE', 0, 'Flag',
               'This variant was used to build the positive training set of good variants')

    culprit, good, tot = {}, {}, 0.0
    anno_texts = ['QD', 'FS', 'BaseQRankSum', 'SOR', 'MQRankSum', 'ReadPosRankSum']

    logger.info("Outputting to %s ..." % opt.output_vcf_file_name)
    OUT = Open(opt.output_vcf_file_name, "wb", isbgz=True) if opt.output_vcf_file_name.endswith(".gz") else \
        open(opt.output_vcf_file_name, "w")

    for k, h in sorted(h_info.header.items(), key=lambda d: d[0]):
        OUT.write("\n".join(h) + "\n")

    n, j, monitor = 0, 0, True
    with Open(opt.vcf_infile, 'r') as I:
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

            # get INFO
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
            OUT.write('\t'.join(col) + "\n")

    OUT.close()

    if opt.output_vcf_file_name.endswith(".gz"):
        # create tabix index
        tabix_index(opt.output_vcf_file_name, force=True, seq_col=0, start_col=1, end_col=1)

    logger.info('Finish Outputting %d lines.\n' % n)

    ## Output Summary
    logger.info('[Summmary] Here is the summary information:')
    for k, v in sorted(good.items(), key=lambda k: k[0]):
        logger.info(('  ** Variant Site score >= %d: %d\t%0.2f' %
                     (k, v, v * 100 / tot)))

    for k, v in sorted(culprit.items(), key=lambda k: k[0]):
        logger.info(('  ** Culprit by %s: %d\t%.2f' % (k, v, v * 100.0 / tot)))


def apply_VQSR(opt):
    """Apply a score cutoff to filter variants."""

    logger.info("Find a VQSLOD cutoff base on %.2f truth set sensitivity level "
                "... ..." % opt.truth_sensitivity_level)

    cdef int total_variant_num = 0
    truth_set_vqlod = []
    false_set_vqlod = []
    with Open(opt.vcf_infile, 'r') as I:
        for line in I:
            if line.startswith('#'):
                continue

            total_variant_num += 1
            col = line.strip().split()

            # get INFO
            vcf_info = {}
            for info in col[7].split(';'):
                cc = info.split('=')
                vcf_info[cc[0]] = cc[-1]

            if 'VQSLOD' not in vcf_info:
                logger.error("Missing VQSLOD, may because you have not run VQSR yet. Abort")
                sys.exit(1)

            if 'POSITIVE_TRAIN_SITE' in vcf_info:
                truth_set_vqlod.append(float(vcf_info['VQSLOD']))

            if 'NEGATIVE_TRAIN_SITE' in vcf_info:
                false_set_vqlod.append(float(vcf_info['VQSLOD']))

    # reverse sorted
    truth_set_vqlod.sort(reverse=True)
    truth_set_num = len(truth_set_vqlod)
    ts_index = int(round(opt.truth_sensitivity_level * truth_set_num)) - 1
    if ts_index < 0:
        ts_index = 0

    vqlod_cutoff = truth_set_vqlod[ts_index]

    false_set_num = len(false_set_vqlod)
    if false_set_num == 0:
        false_set_num = -1

    false_num = 0
    for q in sorted(false_set_vqlod, reverse=True):
        if q >= vqlod_cutoff:
            false_num += 1
        else:
            break

    logger.info("VQLOD cutoff is %s for %s truth set sensitivity level and there are %.2f sensitivity false "
                "variants." % (vqlod_cutoff, opt.truth_sensitivity_level, float(false_num) / false_set_num))

    logger.info("Outputting to %s ..." % opt.output_vcf_file_name)
    OUT = Open(opt.output_vcf_file_name, "wb", isbgz=True) if opt.output_vcf_file_name.endswith(".gz") else \
        open(opt.output_vcf_file_name, "w")

    cdef int pass_variant_num = 0
    with Open(opt.vcf_infile, 'r') as I:
        for line in I:

            if line.startswith('#'):
                OUT.write(line.strip()+"\n")
                continue

            col = line.strip().split()

            qd = re.search(r';?VQSLOD=([^;]+)', col[7])
            vqslod = float(qd.group(1))

            if vqslod >= vqlod_cutoff:
                col[6] = "PASS"
                pass_variant_num += 1

            OUT.write('\t'.join(col) + "\n")

    OUT.close()
    logger.info("There are a total of %d variants, %d of which are PASS base on the VQSLOD "
                "cutoff." % (total_variant_num, pass_variant_num))

    if opt.output_vcf_file_name.endswith(".gz"):
        # create tabix index
        tabix_index(opt.output_vcf_file_name, force=True, seq_col=0, start_col=1, end_col=1)

    return
