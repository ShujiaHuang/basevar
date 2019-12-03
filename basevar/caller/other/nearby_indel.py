"""
This module use for extract nearby indel information for each candidate
varaints.


Author: Shujia Huang
Date: 2017-11-06
"""
import sys
import time

import numpy as np

from basevar.caller.vqsr import vcfutils
from basevar.io.openfile import Open
from basevar.io.BGZF.tabix import TabixFile, tabix_index

class NearbyIndel(object):

    def __init__(self, in_vcf_file, in_cvg_file, output_file, nearby_distance):

        self.in_vcf_file = in_vcf_file
        self.in_cvg_tb = TabixFile(in_cvg_file)
        self.nearby_indel_dis = nearby_distance
        self.output_file_name = output_file

    def _close_input_file(self):
        self.in_cvg_tb.close()

    def _region_indel_sdi(self, chr_id, start, end):
        """
        Calculate the diversity of indel by Shannon's diversity index
        https://zh.wikipedia.org/wiki/%E5%A4%9A%E6%A0%B7%E6%80%A7%E6%8C%87%E6%95%B0
        """
        total_indel_num, indel_type = 0, {}
        for r in self.in_cvg_tb.fetch(chr_id, start=start-1, end=end):
            # chrM    30      T       16150   5       8       2       16135   +1C|1,+AA|2   -0.0    7302,8833,4,4
            col = r.strip().split()
            if col[8] == '.':
                continue

            for indel in [s.upper() for s in col[8].split(',')]:
                indel, n = indel.split('|')
                n = int(n)

                total_indel_num += n
                indel_type[indel] = indel_type.get(indel, 0) + n

        # Indel species
        i_sp = len(indel_type)

        # calculate the SDI
        sdi = np.sum([-1.0*(float(v)/total_indel_num) * np.log2(float(v)/total_indel_num)
                      for k, v in indel_type.items()]) if indel_type else 0.0

        return i_sp, total_indel_num, round(sdi, 3)

    def run(self):

        # Get VCF Header
        h_info = vcfutils.Header()
        with Open(self.in_vcf_file, 'r') as I:
            for line in I:
                # Record the header information
                if line.startswith("#"):
                    h_info.record(line.strip())
                else:
                    break

        h_info.add("INFO", "Indel_SDI", 1, "Float", "Indel diversity by Shannon's diversity index. The less the better.")
        h_info.add("INFO", "Indel_SP", 1, "Integer", "Indel species around this position. The less the better.")
        h_info.add("INFO", "Indel_TOT", 1, "Integer", "Number of Indel around this position. The less the better.")

        # Final output file
        if self.output_file_name == "-":
            OUT = sys.stdout
        else:
            OUT = Open(self.output_file_name, 'wb', isbgz=True if self.output_file_name.endswith(".gz") else False)

        for k, h in sorted(h_info.header.items(), key=lambda d: d[0]):
            OUT.write("\n".join(h) + "\n")

        with Open(self.in_vcf_file, 'r') as I:

            n, monitor = 0, True
            for r in I:

                if r.startswith('#'):
                    continue

                n += 1
                if n % 100000 == 0:
                    sys.stderr.write('** Output lines %d %s\n' % (n, time.asctime()))

                col = r.strip().split()
                pos = int(col[1])  # Chromsome position

                # Deal with the INFO line
                vcfinfo = {}
                for info in col[7].split(';'):
                    k = info.split('=')[0]

                    if monitor and k in vcfinfo:
                        monitor = False
                        sys.stderr.write(('[WARNING] The tag: %s double hits in the INFO column at %s\n' %
                                          (k, self.in_vcf_file)))
                    vcfinfo[k] = info

                start = pos - self.nearby_indel_dis if pos > self.nearby_indel_dis else 1
                end = pos + self.nearby_indel_dis

                indel_sp, indel_tot, indel_sdi = self._region_indel_sdi(col[0], start, end)
                vcfinfo['Indel_SDI'] = 'Indel_SDI=' + str(indel_sdi)
                vcfinfo['Indel_SP'] = 'Indel_SP=' + str(indel_sp)
                vcfinfo['Indel_TOT'] = 'Indel_TOT=' + str(indel_tot)

                col[7] = ';'.join(sorted(vcfinfo.values()))
                OUT.write("%s\n" % "\t".join(col))

        self._close_input_file()
        OUT.close()

        if self.output_file_name.endswith(".gz"):
            # create tabix index
            tabix_index(self.output_file_name, force=True, seq_col=0, start_col=1, end_col=1)

        return self
