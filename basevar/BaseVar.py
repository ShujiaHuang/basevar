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
    basetype_cmd.add_argument('-I', '--input', dest='input', metavar='BAM/CRAM', action='append', default=[],
                              help='BAM/SAM/CRAM file containing reads. This argument could be specified at '
                                   'least once.')
    basetype_cmd.add_argument('-L', '--aligne-file-list', dest='infilelist', metavar='BamfilesList',
                              help='BAM/CRAM files list, one file per row.')
    basetype_cmd.add_argument('-R', '--reference', dest='referencefile', metavar='Reference_fasta', required=True,
                              help='Input reference fasta file.', default='')

    basetype_cmd.add_argument('-q', dest='mapq', metavar='INT', type=int, default=10,
                              help='Only include reads with mapping quality >= INT. [10]')

    basetype_cmd.add_argument('--output-vcf', dest='outvcf', type=str,
                              help='Output VCF file. If not provide will just output position coverage.')
    basetype_cmd.add_argument('--output-cvg', dest='outcvg', type=str, required=True,
                              help='Output position coverage file.')

    basetype_cmd.add_argument('--positions', metavar='position-file', type=str, dest='positions',
                              help='skip unlisted positions (chrid pos). -p and --region could be provided '
                                   'simultaneously.')
    basetype_cmd.add_argument('--regions', metavar='chr:start-end', type=str, dest='regions', default='',
                              help='Skip positions which not in these regions. Comma delimited list of regions '
                                   '(chr:start-end). Could be a file contain the regions. This parameter could '
                                   'be provide with -L simultaneously')

    # The number of output subfiles
    basetype_cmd.add_argument('-B', '--batch-count', dest='batchcount', metavar='NUM', type=int, default=200,
                              help='Simple size per batch. [200]')

    basetype_cmd.add_argument('--nCPU', dest='nCPU', metavar='INT', type=int, default=1,
                              help='Number of processer to use. [1]')
    basetype_cmd.add_argument('-m', '--min-af', dest='min_af', type=float, metavar='float', default=None,
                              help='Setting prior precision of MAF and skip uneffective caller positions. Usually '
                                   'you can set it to be min(0.001, 100/x), x is the number of your input BAM files.'
                                   '[min(0.001, 100/x, cmm.MINAF)]. '
                                   'Probably you don\'t need to take care about this parameter.')

    # special parameter for calculating specific population allele frequence
    basetype_cmd.add_argument('--pop-group', dest='pop_group_file', metavar='GroupListFile', type=str,
                              help='Calculating the allele frequency for specific population.')

    basetype_cmd.add_argument('--filename-has-samplename', dest='filename_has_samplename', action='store_true',
                              help="Sample id should be the first element in filename and been separated by '.' . "
                                   "This will save a lot of time if you have thousands of bamfiles.")

    # smart rerun
    basetype_cmd.add_argument('--smart-rerun', dest='smartrerun', action='store_true',
                              help='Re-run basetype process by checking batchfiles when you set this parameter.')

    # VQSR commands
    vqsr_cmd = commands.add_parser('VQSR', help='Variants Recalibrator')
    vqsr_cmd.add_argument('-I', '--input', dest='vcfInfile', metavar='VCF', action='append', required=True,
                          help='Input VCF file. This argument should be specified at least once.')
    vqsr_cmd.add_argument('-T', '--Train', dest='trainData', metavar='VCF', required=True,
                          help='Traning data set at true category.')
    vqsr_cmd.add_argument('-f', '--fig', dest='figure', metavar='FIG',
                          help='The prefix of figure. [figout]', default='figout')

    # For Coverage
    coverage_cmd = commands.add_parser('coverage', help='Calculating coverage depth for the whole genome '
                                                        'or given regions/positions')
    coverage_cmd.add_argument('-L', '--aligne-file-list', dest='infilelist', metavar='FILE', required=True,
                              help='Input alignmernt file list.', default='')
    coverage_cmd.add_argument('-R', '--reference', dest='referencefile', metavar='FILE', required=True,
                              help='Input reference fasta file.')
    coverage_cmd.add_argument('-O', '--outputfile', dest='outputfile', metavar='FILE',
                              default='out', help='Output file. [out]')

    coverage_cmd.add_argument('-P', '--positions', metavar='position-file', dest='positions',
                              help='skip unlisted positions (chr pos)', default='')
    coverage_cmd.add_argument('--regions', metavar='chr_id:start-end', dest='regions',
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
    }

    args = parser_commandline_args()

    sys.stderr.write('\n** %s Start at %s **\n\n' % (args.command, time.asctime()))

    runner[args.command](args)

    elasped_time = datetime.now() - START_TIME
    sys.stderr.write('** %s done at %s, %d seconds elapsed **\n' % (
        args.command, time.asctime(), elasped_time.seconds))


if __name__ == '__main__':
    main()
