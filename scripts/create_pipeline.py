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
    exe_prog = pardir + '/basevar/BaseVar.py basetype'

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
                      help='Set specific region size', default=100000)
    optp.add_argument('-c', '--chrom', metavar='STR', dest='chrom',
                      help='skip comma delimited unlisted chrom. e.g chr1,chr2',
                      default='')
    optp.add_argument('--nCPU', dest='nCPU', metavar='INT', type=int,
                      help='Number of processer to use. [1]', default=1)

    optp.add_argument('-b', '--bgzip', metavar='STR', dest='bgzip',
                      help='The path of bgzip',
                      default='/home/huangshujia/local/bin/bgzip')
    optp.add_argument('-t', '--tabix', metavar='STR', dest='tabix',
                      help='The path of tabix',
                      default='/home/huangshujia/local/bin/tabix')

    ## Parameters for BaseType
    optp.add_argument('-L', '--align-file-list', dest='infilelist', metavar='FILE',
                      help='BAM/CRAM files list, one file per row.', default='')
    optp.add_argument('-m', '--min_af', dest='min_af', type=float, metavar='MINAF', default=0.001,
                      help='By setting min AF to skip uneffective caller positions to accelerate program speed. [0.001]')

    opt = optp.parse_args()
    opt.delta = int(opt.delta)
    opt.infilelist = os.path.abspath(opt.infilelist)
    opt.outdir = os.path.abspath(opt.outdir)

    bgzip = opt.bgzip
    tabix = opt.tabix

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
            print(f'time python {exe_prog} --nCPU {opt.nCPU} -L {opt.infilelist} --regions {reg} '
                  f'-m {opt.min_af} --reference {opt.reference} '
                  f'--output-vcf {opt.outdir}/{outfile_prefix}.vcf '
                  f'--output-cvg {opt.outdir}/{outfile_prefix}.cvg.tsv --smart-rerun && '
                  f'{bgzip} -f {opt.outdir}/{outfile_prefix}.vcf && '
                  f'{bgzip} -f {opt.outdir}/{outfile_prefix}.cvg.tsv && '
                  f'{tabix} -f -p vcf {opt.outdir}/{outfile_prefix}.vcf.gz && '
                  f'{tabix} -f -b 2 -e 2 {opt.outdir}/{outfile_prefix}.cvg.tsv.gz && '
                  f'echo "** {outfile_prefix} done **"')


if __name__ == '__main__':
    creat_basetype_pipe()

