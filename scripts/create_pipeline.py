"""
Author: Shujia Huang
Date: 2017-01-19
"""
import os
import sys
import argparse

def load_reference_fai(in_fai, chroms=None):

    ref = []
    with open(in_fai) as fh:

        for r in fh:
            # chr2  243199373  254235646  50 51
            col = r.strip().split()
            if chroms != None and len(chroms):
                if col[0] in chroms:
                    ref.append([col[0], 1, int(col[1])])
            else:
                ref.append([col[0], 1, int(col[1])])

    return ref


def creat_basetype_pipe():

    pardir = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), os.path.pardir))
    exe_prog = pardir + '/bin/basevar basetype'

    optp = argparse.ArgumentParser()
    optp.add_argument('-o', '--outdir', metavar='STR', dest='outdir',
                      help='The output directory', default='')
    optp.add_argument('-R', '--reference_fasta', metavar='FILE', dest='reference',
                      help='The reference fa file', default='')
    optp.add_argument('-f', '--ref_fai', metavar='FILE', dest='ref_fai',
                      help='The reference fai file', default='')
    optp.add_argument('-r', '--regions', metavar='Region', dest='regions',
                      help='skip positions not in (chr:start-end)', default='')
    optp.add_argument('-d', '--delta', metavar='INT', dest='delta',
                      help='Set specific region size', default=2000000)
    optp.add_argument('-c', '--chrom', metavar='STR', dest='chrom',
                      help='skip comma delimited unlisted chrom. e.g chr1,chr2',
                      default='')
    optp.add_argument('-t', '--thread', dest='n_thread', metavar='INT', type=int,
                      help='Number of processer to use. [1]', default=1)

    ## Parameters for BaseType
    optp.add_argument('-L', '--align-file-list', dest='infilelist', metavar='FILE',
                      help='BAM/CRAM files list, one file per row.', default='')
    optp.add_argument('-m', '--min_af', dest='min_af', type=float, metavar='MINAF', default=0.001,
                      help='By setting min AF to skip uneffective caller positions to accelerate program speed. [0.001]')
    optp.add_argument('-G', '--pop-group', dest='pop_group', metavar='FILE',
                      help='Calculating the allele frequency for specific population.', default='')

    opt = optp.parse_args()
    opt.delta = int(opt.delta)
    opt.infilelist = os.path.abspath(opt.infilelist)
    opt.outdir = os.path.abspath(opt.outdir)

    chroms = opt.chrom.strip().split(',') if opt.chrom else []
    ref_fai = load_reference_fai(opt.ref_fai, chroms) if opt.ref_fai else []

    if len(opt.regions):
        # Todo: 可以增加多 regions 的功能，逗号分开，做数组即可
        chrid, reg = opt.regions.strip().split(':')
        reg = list(map(int, reg.split('-')))
        ref_fai[chrid, reg[0], reg[1]] if (not chroms) or (chrid in chroms) else []

    for chr_id, reg_start, reg_end in ref_fai:
        for i in range(reg_start-1, reg_end, opt.delta):
            start = i + 1
            end = i + opt.delta if i + opt.delta <= reg_end else reg_end
            reg = chr_id + ':' + str(start) + '-' + str(end)

            outfile_prefix = chr_id + '_' + str(start) + '_' + str(end)
            if opt.pop_group:
                print(f'time {exe_prog} -t {opt.n_thread} -R {opt.reference} -L {opt.infilelist} '
                      f'-G {opt.pop_group} -r {reg} --min-af={opt.min_af} '
                      f'--output-vcf {opt.outdir}/{outfile_prefix}.vcf.gz '
                      f'--output-cvg {opt.outdir}/{outfile_prefix}.cvg.tsv.gz '
                      f'--smart-rerun > {opt.outdir}/{outfile_prefix}.log && '
                      f'echo "** {outfile_prefix} done **"')
            else:    
                print(f'time {exe_prog} -t {opt.n_thread} -R {opt.reference} -L {opt.infilelist} '
                      f'-r {reg} --min-af={opt.min_af} '
                      f'--output-vcf {opt.outdir}/{outfile_prefix}.vcf.gz '
                      f'--output-cvg {opt.outdir}/{outfile_prefix}.cvg.tsv.gz '
                      f'--smart-rerun > {opt.outdir}/{outfile_prefix}.log && '
                      f'echo "** {outfile_prefix} done **"')


if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print("Please type: -h or --help to show the help message.\n", file=sys.stderr)
        sys.exit(1)

    creat_basetype_pipe()

