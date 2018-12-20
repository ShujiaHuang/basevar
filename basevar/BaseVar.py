"""
This is the main program of BaseVar. It's the toppest layer of
BaseVar's tool sets.

Autor: Shujia Huang
Date: 2016-10-06 16:38:00
"""
import argparse
import sys
import time

from datetime import datetime


def basetype(args):
    from caller.executor import BaseTypeBamRunner
    bt = BaseTypeBamRunner(args)
    bt.run()

    return


# def vqsr():
#     from caller.executor import VQSRRuner
#     vq = VQSRRuner()
#     vq.run()
#
#     return


def nearby_indel(args):
    from caller.executor import NearbyIndelRunner
    nbi = NearbyIndelRunner(args)
    nbi.run()

    return


def merge(args):
    from caller.executor import MergeRunner

    mg = MergeRunner(args)
    mg.run()

    return


def coverage(args):
    from caller.executor import CoverageRunner
    cvg = CoverageRunner(args)
    cvg.run()

    return


def fusion():
    from caller.executor import FusionRunner
    cf = FusionRunner()
    cf.run()

    return


def fusionbasetype():
    from caller.executor import BaseTypeFusionRunner
    ft = BaseTypeFusionRunner()
    ft.run()

    return


def parser_commandline_args():
    desc = "BaseVar: A python software for calling variants without calling genotype."
    cmdparse = argparse.ArgumentParser(description=desc)
    commands = cmdparse.add_subparsers(dest="command", title="BaseVar Commands")

    basetype_cmd = commands.add_parser('basetype', help='Variants Caller')
    basetype_cmd.add_argument('-I', '--aligne-file-list', dest='infilelist', metavar='Bamfiles', required=True,
                              help='BAM/CRAM file list, one line per file.')
    basetype_cmd.add_argument('-R', '--reference', dest='referencefile', metavar='Reference_fasta', required=True,
                              help='Input reference fasta file.', default='')
    basetype_cmd.add_argument('-O', '--outprefix', dest='outprefix', metavar='VCF_Prefix',
                              default='out', help='The prefix of output files. [out]')

    basetype_cmd.add_argument('-L', '--positions', metavar='positions', type=str, dest='positions', default='',
                              help='skip unlisted positions (chrid pos). -L and --region could be provided '
                                   'simultaneously. [None]')
    basetype_cmd.add_argument('--region', metavar='chr:start-end', type=str, dest='regions', default='',
                              help='Skip positions which not in these regions. Comma delimited list of regions '
                                   '(chr:start-end). Could be a file contain the regions. This parameter could '
                                   'be provide with -L simultaneously')

    basetype_cmd.add_argument('-q', dest='mapq', metavar='INT', type=int, default=10,
                              help='Only include reads with mapping quality >= INT. [10]')
    # The number of output subfiles
    basetype_cmd.add_argument('--batch-count', dest='batchcount', metavar='NUM', type=int, default=1000,
                              help='Number of samples in a batch file. [1000]')

    basetype_cmd.add_argument('--nCPU', dest='nCPU', metavar='Int', type=int, default=1,
                              help='Number of processer to use. [1]')
    basetype_cmd.add_argument('-m', '--min_af', dest='min_af', type=float, metavar='float', default=None,
                              help='Setting prior precision of MAF and skip uneffective caller positions. Usually '
                                   'you can set it to be min(0.001, 100/x), x is the number of your input BAM files.'
                                   '[min(0.001, 100/x, cmm.MINAF)]. Probably you donot need to care about this parameter.')

    # special parameter for calculating specific population allele frequence
    basetype_cmd.add_argument('--pop-group', dest='pop_group_file', metavar='Group-List-File', type=str,
                              help='Calculating the allele frequency for specific population.')

    basetype_cmd.add_argument('--filename-has-samplename', dest='filename_has_samplename', type=bool, default=False,
                              help="Sample id should be the first element in filename and been separated by '.' . "
                                   "This will save a lot of time if you have thousands of bamfiles. [False]")

    # smart re-run
    basetype_cmd.add_argument('--smart-rerun', dest='smartrerun', action='store_true',
                              help='Re-run basetype process by checking batchfiles '
                                   'when you set this parameter.')

    # special parameter to limit the function of BaseType
    basetype_cmd.add_argument('--justdepth', dest='justdepth', action='store_true',
                              help='Setting for just output genome coverage.')

    # VQSR commands
    vqsr_cmd = commands.add_parser('VQSR', help='Variants Recalibrator')
    vqsr_cmd.add_argument('-i', '--InVcf', dest='vcfInfile', metavar='VCF',
                          help='VCF for predict.', default=[])
    vqsr_cmd.add_argument('-T', '--Train', dest='trainData', metavar='TRU',
                          help='Traing data set at true category', default=[])
    vqsr_cmd.add_argument('-f', '--fig', dest='figure', metavar='FIG',
                          help='The prefix of figure.', default='figtest')

    # For Coverage
    coverage_cmd = commands.add_parser('coverage', help='Calculating coverage depth for the whole genome '
                                                        'or given regions/positions')
    coverage_cmd.add_argument('-l', '--aligne-file-list', dest='infilelist', metavar='FILE', required=True,
                              help='Input alignmernt file list.', default='')
    coverage_cmd.add_argument('-r', '--reference', dest='referencefile', metavar='FILE', required=True,
                              help='Input reference fasta file.')
    coverage_cmd.add_argument('-O', '--outputfile', dest='outputfile', metavar='FILE',
                              default='out', help='Output file. [out]')

    coverage_cmd.add_argument('-P', '--positions', metavar='FILE', dest='positions',
                              help='skip unlisted positions (chr pos)', default='')
    coverage_cmd.add_argument('-R', '--regions', metavar='chr:start-end', dest='regions',
                              help='skip positions not in (chr:start-end)', default='')

    coverage_cmd.add_argument('--nCPU', dest='nCPU', metavar='INT', type=int,
                              help='Number of processer to use. [1]', default=1)

    # Merge files
    merge_cmd = commands.add_parser('merge', help='Merge bed/vcf files')
    merge_cmd.add_argument('-L', '--file-list', dest='infilelist', metavar='FILE', required=True,
                           help='Input files\' list.', default='')
    merge_cmd.add_argument('-O', '--outputfile', dest='outputfile', metavar='FILE', default='out',
                           help='Output file. [out]')

    # Add nearby indels for variants
    nbi_cmd = commands.add_parser('NearByIndel', help='Calculating and adding Nearby Indel density and '
                                                      'indel type information for each variants in VCF')
    nbi_cmd.add_argument('-i', '--in-vcf-file', dest='in_vcf_file', metavar='VCF_FILE', required=True,
                         help='The input vcf files')
    nbi_cmd.add_argument('-c', '--in-cvg-file', dest='in_cvg_file', metavar='BaseVar_CVG_FILE', required=True,
                         help='Input coverage file which has indel information')
    nbi_cmd.add_argument('-d', '--nearby-distance-around-indel', dest='nearby_dis_around_indel', type=int,
                         help='The distance around indel. [16]', default=16)

    return cmdparse.parse_args()


def main():
    START_TIME = datetime.now()

    runner = {
        'basetype': basetype,
        'merge': merge,
        'coverage': coverage,
        'NearByIndel': nearby_indel,
        # 'VQSR': vqsr
        # 'fusion': fusion,
        # 'fusionbasetype': fusionbasetype,
    }

    args = parser_commandline_args()

    sys.stderr.write('** %s Start at %s **\n\n' % (args.command, time.asctime()))

    runner[args.command](args)

    elasped_time = datetime.now() - START_TIME
    sys.stderr.write('** %s done at %s, %d seconds elapsed **\n' % (
        args.command, time.asctime(), elasped_time.seconds))


if __name__ == '__main__':
    main()
