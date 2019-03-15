"""
This is the main program of BaseVar. It's the toppest layer of
BaseVar's tool sets.

Autor: Shujia Huang
Date: 2016-10-06 16:38:00
"""
import argparse
import sys
import time


def basetype(args):
    from caller.executor import BaseTypeRunner

    if args.outbatchfile and (args.outvcf or args.outcvg):
        sys.stderr.write("[ERROR] Don't set '--output-vcf' or '--output-cvg' if you have set "
                         "'--output-batch-file'.\n\n")
        sys.exit(1)

    # Make sure you have set at least one output parameter.
    if not (args.outbatchfile or args.outvcf or args.outcvg):
        sys.stderr.write("[ERROR] Missing '--output-vcf' or '--output-cvg' or "
                         "'--output-batch-file'.\n\n")
        sys.exit(1)

    if not args.outbatchfile and not args.outcvg:
        sys.stderr.write("[ERROR] argument '--output-cvg' is required if not set '--output-batch-file'\n\n")
        sys.exit(1)

    if args.smartrerun:
        sys.stderr.write("************************************************\n"
                         "******************* WARNING ********************\n"
                         "************************************************\n"
                         ">>>>>>>> You have setted `smart rerun` <<<<<<<<<\n"
                         "Please make sure all the parameters are the same\n"
                         "with your previous commands.\n"
                         "************************************************\n\n")

    # Make sure you have set at least one bamfile.
    if not args.input and not args.infilelist:
        sys.stderr.write("[ERROR] Missing input BAM/CRAM files.\n\n")
        sys.exit(1)

    # The main function
    bt = BaseTypeRunner(args)

    if not args.outbatchfile:
        processer = bt.basevar_caller()
    else:
        processer = bt.batch_generator()

    is_success = True
    for p in processer:
        if p.exitcode != 0:
            is_success = False

    return is_success


def basetypebatch(args):
    from caller.executor import BaseTypeBatchRunner

    # Make sure you have set at least one input file.
    if not args.input and not args.infilelist:
        sys.stderr.write("[ERROR] Missing input batch files.\n\n")
        sys.exit(1)

    bt = BaseTypeBatchRunner(args)

    processer = bt.basevar_caller()

    is_success = True
    for p in processer:
        if p.exitcode != 0:
            is_success = False

    return is_success


def vqsr(args):
    # Todo: VQSR need to improve
    from caller.executor import VQSRRuner
    vq = VQSRRuner(args)
    vq.run()

    return True


def nearby_indel(args):
    from caller.executor import NearbyIndelRunner
    nbi = NearbyIndelRunner(args)
    nbi.run()

    return True


def merge(args):
    from caller.executor import MergeRunner

    mg = MergeRunner(args)
    mg.run()

    return True


def coverage(args):
    from caller.executor import CoverageRunner
    cvg = CoverageRunner(args)
    processer = cvg.run()

    is_success = True
    for p in processer:
        if p.exitcode != 0:
            is_success = False

    return is_success


def parser_commandline_args():
    desc = "BaseVar: A python software to call variants without calling genotype for ultra low pass " \
           "sequencing data."

    cmdparse = argparse.ArgumentParser(description=desc)
    commands = cmdparse.add_subparsers(dest="command", title="BaseVar Commands")

    basetype_cmd = commands.add_parser('basetype', help='Variants Caller')
    basetype_cmd.add_argument('-I', '--input', dest='input', metavar='BAM/CRAM', action='append', default=[],
                              help='BAM/SAM/CRAM file containing reads. This argument could be specified at '
                                   'least once.')
    basetype_cmd.add_argument('-L', '--align-file-list', dest='infilelist', metavar='BamfilesList',
                              help='BAM/CRAM files list, one file per row.')
    basetype_cmd.add_argument('-R', '--reference', dest='referencefile', metavar='Reference_fasta', required=True,
                              help='Input reference fasta file.')

    basetype_cmd.add_argument('-q', dest='mapq', metavar='INT', type=int, default=10,
                              help='Only include reads with mapping quality >= INT. [10]')

    basetype_cmd.add_argument('--output-vcf', dest='outvcf', type=str,
                              help='Output VCF file. If not provide will skip variants discovery and just output '
                                   'position coverage file which filename is provided by --output-cvg.')
    basetype_cmd.add_argument('--output-cvg', dest='outcvg', type=str, help='Output position coverage file.')
    basetype_cmd.add_argument('--output-batch-file', dest='outbatchfile', type=str, help='Just emitting batch-file.')

    basetype_cmd.add_argument('--positions', metavar='position-list-file', type=str, dest='positions',
                              help='skip unlisted positions one per row. The position format in the file could '
                                   'be (chrid pos) and (chrid start end) in mix. This parameter could be used '
                                   'with --regions simultaneously.')
    basetype_cmd.add_argument('--regions', metavar='chr:start-end', type=str, dest='regions', default='',
                              help='Skip positions which not in these regions. This parameter could be a list of '
                                   'comma deleimited genome regions(e.g.: chr:start-end,chr:start-end) or a file '
                                   'contain the list of regions. This parameter could be used with --positions '
                                   'simultaneously')

    # The number of output subfiles
    basetype_cmd.add_argument('-B', '--batch-count', dest='batchcount', metavar='INT', type=int, default=200,
                              help='INT simples per batchfile. [200]')
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
                              help="If the name of bamfile is something like 'SampleID.xxxx.bam', "
                                   "you can set this parameter to save a lot of time during get the "
                                   "sample id from BAM header.")

    basetype_cmd.add_argument('--smart-rerun', dest='smartrerun', action='store_true',
                              help='Rerun process by checking batchfiles.')

    # For discovery variants from batchfiles
    btb_cmd = commands.add_parser('basetypebatch',
                                  help='Variants discovery on one or more samples pre-call by basetype')
    btb_cmd.add_argument('-I', '--input', dest='input', metavar='BatchFile', action='append', default=[],
                         help='Input batchfile pre-call by basetype and must be compressed by bgzip algorithm. '
                              'This argument could be specified at least once.')
    btb_cmd.add_argument('-L', '--batch-file-list', dest='infilelist', metavar='BatchfilesList',
                         help='batchfiles list pre-call by basetype and must be compressed by bgzip algorithm. '
                              'One file per row.')
    btb_cmd.add_argument('-R', '--reference', dest='referencefile', metavar='Reference_fasta', required=True,
                         help='Input reference fasta file.')

    btb_cmd.add_argument('--output-vcf', dest='outvcf', type=str,
                         help='Output VCF file. If not provide will skip variants discovery and just output '
                              'position coverage file which filename is provided by --output-cvg.')
    btb_cmd.add_argument('--output-cvg', dest='outcvg', type=str, required=True,
                         help='Output position coverage file.')

    btb_cmd.add_argument('--positions', metavar='position-list-file', type=str, dest='positions',
                         help='Skip unlisted positions one per row. The position format in the file could '
                              'be (chrid pos) and (chrid start end) in mix. This parameter could be used '
                              'with --regions simultaneously.')
    btb_cmd.add_argument('--regions', metavar='chr:start-end', type=str, dest='regions', default='',
                         help='Skip positions which not in these regions. This parameter could be a list of '
                              'comma deleimited genome regions(e.g.: chr:start-end,chr:start-end) or a file '
                              'contain the list of regions. This parameter could be used with --positions '
                              'simultaneously')

    btb_cmd.add_argument('-m', '--min-af', dest='min_af', type=float, metavar='float', default=None,
                         help='Setting prior precision of MAF and skip uneffective caller positions. Usually '
                              'you can set it to be min(0.001, 100/x), x is the number of your input BAM files.'
                              '[min(0.001, 100/x, cmm.MINAF)]. Probably you don\'t need to take care about this '
                              'parameter.')

    # special parameter for calculating specific population allele frequence
    btb_cmd.add_argument('--pop-group', dest='pop_group_file', metavar='GroupListFile', type=str,
                         help='Calculating the allele frequency for specific population.')
    btb_cmd.add_argument('--nCPU', dest='nCPU', metavar='INT', type=int, default=1,
                         help='Number of processer to use. [1]')

    # VQSR commands
    vqsr_cmd = commands.add_parser('VQSR', help='Variants Recalibrator')
    vqsr_cmd.add_argument('-I', '--input', dest='vcfInfile', metavar='VCF', required=True,
                          help='Input VCF file. This argument should be specified at least once.')
    vqsr_cmd.add_argument('-T', '--Train', dest='trainData', metavar='VCF', required=True,
                          help='Traning data set at true category.')
    vqsr_cmd.add_argument('--fig', dest='figure', metavar='FIG', required=True,
                          help='The prefix of figure.')

    # For Coverage
    coverage_cmd = commands.add_parser('coverage', help='Calculating coverage depth for the whole genome '
                                                        'or given regions/positions')
    coverage_cmd.add_argument('-L', '--align-file-list', dest='infilelist', metavar='FILE', required=True,
                              help='Input alignmernt file list.', default='')
    coverage_cmd.add_argument('-R', '--reference', dest='referencefile', metavar='FILE', required=True,
                              help='Input reference fasta file.')
    coverage_cmd.add_argument('-O', '--outputfile', dest='outputfile', metavar='FILE', default='out',
                              help='Output file. [out]')

    coverage_cmd.add_argument('--positions', metavar='position-list-file', type=str, dest='positions',
                              help='skip unlisted positions one per row. The position format in the file could '
                                   'be (chrid pos) and (chrid  start  end) in mix. This parameter could be used '
                                   'with --regions simultaneously.')
    coverage_cmd.add_argument('--regions', metavar='chr:start-end', type=str, dest='regions', default='',
                              help='Skip positions which not in these regions. This parameter could be a list of '
                                   'comma deleimited genome regions(e.g.: chr:start-end,chr:start-end) or a file '
                                   'contain the list of regions. This parameter could be used with --positions '
                                   'simultaneously')

    coverage_cmd.add_argument('--nCPU', dest='nCPU', metavar='INT', type=int,
                              help='Number of processor to use. [1]', default=1)

    # Merge files
    merge_cmd = commands.add_parser('merge', help='Merge bed/vcf files')
    merge_cmd.add_argument('-I', '--input', dest='input', metavar='FILE', action='append', default=[],
                           help='Input files. This argument could be specify at least once')
    merge_cmd.add_argument('-L', '--file-list', dest='infilelist', metavar='FILE', help='Input files\' list.')
    merge_cmd.add_argument('-O', '--outputfile', dest='outputfile', metavar='FILE', required=True,
                           help='Output file')

    # Add nearby indels for variants
    nbi_cmd = commands.add_parser('NearByIndel', help='Calculating and adding Nearby Indel density and '
                                                      'indel type information for each variants in VCF.')
    nbi_cmd.add_argument('-I', '--in-vcf-file', dest='in_vcf_file', metavar='VCF_FILE', required=True,
                         help='The input vcf files.')
    nbi_cmd.add_argument('-C', '--in-cvg-file', dest='in_cvg_file', metavar='BaseVar_CVG_FILE', required=True,
                         help='Input coverage file which has indel information.')
    nbi_cmd.add_argument('-D', '--nearby-distance-around-indel', dest='nearby_dis_around_indel', metavar='INT',
                         type=int, default=16, help='The distance around indels. [16]')

    return cmdparse.parse_args()


def main():
    start_time = time.time()
    runner = {
        'basetype': basetype,
        'basetypebatch': basetypebatch,
        'merge': merge,
        'coverage': coverage,
        'NearByIndel': nearby_indel,
        'VQSR': vqsr
    }

    args = parser_commandline_args()
    sys.stderr.write('\n** %s Start at %s **\n' % (args.command, time.asctime()))

    is_success = runner[args.command](args)

    elapsed_time = time.time() - start_time
    if is_success:
        sys.stderr.write('** %s done at %s, %d seconds elapsed **\n' % (
            args.command, time.asctime(), elapsed_time))
    else:
        sys.stderr.write('[ERROR] Catch some exception on %s, so "%s" is not done, %d seconds elapsed\n' % (
            time.asctime(), args.command, elapsed_time))
        sys.exit(1)


if __name__ == '__main__':
    main()
