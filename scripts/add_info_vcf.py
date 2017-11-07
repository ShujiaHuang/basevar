"""
Author: Shujia Huang
Date: 2017-11-07
"""
import sys
import time
import gzip


def main(argv):

    target_file, in_vcf_file = argv

    info = {}
    with open(target_file) as I:

        # scan the file
        for r in I:
            # chr   pos   value

            if r.startswith('#'): continue

            col = r.strip().split()
            info[col[0]+':'+col[1]] = col[2]

    with gzip.open(in_vcf_file) as I:

        n, monitor = 0, True
        for r in I:
            if r.startswith('##FORMAT=<ID=GT,'):
                print ('##INFO=<ID=DM,Number=1,Type=Float,'
                       'Description="Differencial mapping index">')

            if r.startswith('#'):
                print (r.strip())
                continue

            n += 1
            if n % 100000 == 0:
                sys.stderr.write('** Output lines %d %s\n' %
                                 (n, time.asctime()))

            col = r.strip().split()

            # Deal with the INFO line
            vcfinfo = {}
            for info in col[7].split(';'):
                k = info.split('=')[0]

                if monitor and k in vcfinfo:
                    monitor = False
                    sys.stderr.write(('[WARNING] The tag: %s double hits in '
                                      'the INFO column at %s\n' %
                                      (k, in_vcf_file)))
                vcfinfo[k] = info

            pos_key = col[0] + ':' + col[1]
            vcfinfo['DM'] = 'DM=' + str(info[pos_key])

            col[7] = ';'.join(sorted(vcfinfo.values()))
            print ('\t'.join(col))


if __name__ == '__main__':

    main(sys.argv[1:])





