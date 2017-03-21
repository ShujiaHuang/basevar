"""
Author: Shujia Huang
Date: 2017-01-19
"""
import os
import sys
import argparse

def load_reference_fai(in_fai, chroms=None):

    ref = {}
    with open(in_fai) as fh:

        for r in fh:
            # chr2  243199373  254235646  50 51
            col = r.strip().split()
            if chroms != None and len(chroms):
                if col[0] in chroms:
                    ref[col[0]] = [1, int(col[1])]
            else:
                ref[col[0]] = [1, int(col[1])]

    return ref


def creat_basetype_pipe():

    pardir = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), os.path.pardir))
    exe_prog = pardir + '/basevar/BaseVar.py basetype'

    optp = argparse.ArgumentParser()
    optp.add_argument('-f', '--ref_fai', metavar='FILE', dest='ref_fai',
                      help='The reference fai file', default='')
    optp.add_argument('-r', '--regions', metavar='Region', dest='regions',
                      help='skip positions not in (chr:start-end)', default='')
    optp.add_argument('-d', '--delta', metavar='INT', dest='delta',
                      help='Set the region size', default=10000)
    optp.add_argument('-c', '--chrom', metavar='STR', dest='chrom',
                      help='skip comma delimited unlisted chrom. e.g chr1,chr2',
                      default='')
    optp.add_argument('-o', '--outdir', metavar='STR', dest='outdir',
                      help='The output directory', default='')

    optp.add_argument('-b', '--bgzip', metavar='STR', dest='bgzip',
                      help='The path of bgzip',
                      default='/home/huangshujia/local/bin/bgzip')
    optp.add_argument('-t', '--tabix', metavar='STR', dest='tabix',
                      help='The path of tabix',
                      default='/home/huangshujia/local/bin/tabix')

    ## Parameters for BaseType
    optp.add_argument('-l', '--mpileup-list', dest='infilelist', metavar='FILE',
                      help='The input mpileup file list.', default='')
    optp.add_argument('-s', '--sample-list', dest='samplelistfile',
                      metavar='FILE', help='The sample list.')

    opt = optp.parse_args()
    opt.delta = int(opt.delta)
    opt.infilelist = os.path.abspath(opt.infilelist)
    opt.samplelistfile = os.path.abspath(opt.samplelistfile)
    opt.outdir = os.path.abspath(opt.outdir)

    bgzip = opt.bgzip
    tabix = opt.tabix

    chroms = opt.chrom.strip().split(',') if opt.chrom else []
    ref_fai = load_reference_fai(opt.ref_fai, chroms) if opt.ref_fai else {}

    if len(opt.regions):
        chrid, reg = opt.regions.strip().split(':')
        reg = map(int, reg.split('-'))
        ref_fai[chrid] = [reg[0], reg[1]] if (not chroms) or (chrid in chroms) else {}

    for chr_id, (reg_start, reg_end) in sorted(ref_fai.items(), key=lambda x:x[0]):
        for i in range(reg_start-1, reg_end, opt.delta):
            start = i + 1
            end = i + opt.delta if i + opt.delta <= reg_end else reg_end
            reg = chr_id + ':' + str(start) + '-' + str(end)

            outfile_prefix = chr_id + '_' + str(start) + '_' + str(end)
            print ' '.join(['time python '+ exe_prog,
                            '-R '+ reg,
                            '-l '+ opt.infilelist,
                            '-s '+ opt.samplelistfile,
                            '-o '+ opt.outdir + '/' + outfile_prefix,
                            '&& '+ bgzip + ' -f' + opt.outdir + '/' + outfile_prefix + '.vcf',
                            '&& '+ bgzip + ' -f' + opt.outdir + '/' + outfile_prefix + '.cvg.tsv',
                            '&& '+ tabix + ' -f -p vcf' + opt.outdir + '/' + outfile_prefix + '.vcf.gz',
                            '&& '+ tabix + ' -f -b 2 -e 2' + opt.outdir + '/' + outfile_prefix + '.cvg.tsv.gz',
                            '&& echo "** %s done **"' % outfile_prefix])



if __name__ == '__main__':

    creat_basetype_pipe()

