"""
This module contain functions of EM algorithm and Base genotype.

Author: Shujia Huang
Date : 2016-12-16
Update : 2017-01-03

"""
import sys
import multiprocessing

import itertools   # Use the combinations function
import numpy as np
from scipy.stats import chisqprob as sp_chisqprob

import pysam

from . import utils
from . import mpileup
from .algorithm import EM, strand_bias


class BaseType(object):

    def __init__(self, ref_base, bases, quals, cmm=None):
        """A class for calculate the base probability

        Parameters
        ----------

        ``ref_base``: A char, required
            The reference base

        ``bases``: A array like, required
            A list of base for samples.

        ``quals``: An array like, required
            Base quality for ``bases``. The same size with ``bases``
            Cause: The ``quals`` is an integer array which has be converted
                by phred-scale

        ``cmm``: A ``CommonParameter`` class value, required

        """

        self._alt_bases = []
        self._var_qual = 0 # init the variant quality
        self._ref_base = ref_base
        self.cmm = cmm

        #The prior probability for echo base
        self.prior_prob = []
        self.eaf = {} ## estimated allele frequency by EM
        self.depth = {b:0 for b in self.cmm.BASE}

        quals = np.array(quals)
        self.qual_pvalue = 1.0 - np.exp(self.cmm.MLN10TO10 * quals)
        for i, b in enumerate(bases):
            ## Individual per row for [A, C, G, T] and give a prior probability
            if b != 'N':  # ignore all the 'N' base sample
                self.prior_prob.append([self.qual_pvalue[i]
                                        if b == t else (1.0-self.qual_pvalue[i])/3
                                        for t in self.cmm.BASE])
                # record depth for [ACGT]
                if b in self.depth: # ignore '*'
                    self.depth[b] += 1

        self.prior_prob = np.array(self.prior_prob)
        self.total_depth = float(sum(self.depth.values()))

        return

    def ref_base(self):
        return self._ref_base

    def alt_bases(self):
        return self._alt_bases

    def var_qual(self):
        return self._var_qual

    def debug(self):
        print(self.ref_base(), self.alt_bases(),
              self.var_qual(), self.depth, self.eaf)

    def _set_base_likelihood(self, bases):
        """
        init the base likelihood by bases

        ``bases``: a list like
        """
        total_depth = float(sum([self.depth[b] for b in bases]))
        base_likelihood = np.zeros(len(self.cmm.BASE))  # [A, C, G, T] set to 0.0

        if total_depth > 0:
            for b in bases:
                base_likelihood[self.cmm.BASE2IDX[b]] = self.depth[b]/total_depth

        return np.array(base_likelihood)

    def _f(self, bases, n):
        """
        Calculate population likelihood for all the combination of bases

        Parameters
        ----------

        ``bases``: 1d array like
            A list of bases from [A, C, G, T]

        ``n``: Integer
            The combination number. n must less or equal
            to the length of ``bases``

        Return
        ------
        ``bc``: array=like, combination bases
        ``lr``: Likelihood of ``bc``

        Example
        -------

        >>> import itertools
        >>> bases = ['A', 'C', 'G', 'T']
        >>> bc=[i for i in itertools.combinations(bases,3)]
        >>> bc
        ... [('A', 'C', 'G'), ('A', 'C', 'T'), ('A', 'G', 'T'), ('C', 'G', 'T')]

        """
        bc = []
        lr = []
        bp = []

        for b in [i for i in itertools.combinations(bases, n)]:

            init_likelihood = self._set_base_likelihood(b)
            if sum(init_likelihood) == 0: continue  ## No covert

            _, pop_likelihood, base_expected_prob = EM(
                self.prior_prob,
                np.tile(init_likelihood, (self.prior_prob.shape[0],1)))

            bc.append(b)
            lr.append(np.log(pop_likelihood).sum()) # sum the marginal prob
            bp.append(base_expected_prob)

        return bc, lr, bp

    def lrt(self):
        """The main function.
        likelihood ratio test.
        """
        if self.total_depth == 0: return

        # get effective bases which count frequence > self.cmm.MINAF
        bases = [b for b in self.cmm.BASE
                 if self.depth[b]/self.total_depth > self.cmm.MINAF]

        # init
        _, lr_null, bp = self._f(bases, len(bases))

        base_frq = bp[0]
        lr_alt = lr_null[0]

        chi_sqrt_value = 0
        for n in range(1, len(bases))[::-1]:

            bc, lr_null, bp = self._f(bases, n)
            lrt_chivalue = 2.0 * (lr_alt - lr_null)
            i_argmin = np.argmin(lrt_chivalue)

            lr_alt = lr_null[i_argmin]
            chi_sqrt_value = lrt_chivalue[i_argmin]
            if chi_sqrt_value < self.cmm.LRT_THRESHOLD:
                # Take the null hypothesis and continue
                bases = bc[i_argmin]
                base_frq = bp[i_argmin]
            else:
                # Take the alternate hypothesis
                break

        self._alt_bases = [b for b in bases if b != self._ref_base]
        self.eaf = {b : '%f' % round(base_frq[self.cmm.BASE2IDX[b]], 6)
                    for b in bases if b != self._ref_base}

        # calculate the variant quality
        # Todo: improve the calculation method for var_qual
        if len(self._alt_bases):

            if len(bases) == 1 and self.depth[bases[0]] / self.total_depth > 0.5:
                # mono-allelelic
                self._var_qual = 5000

            else:
                chi_prob = sp_chisqprob(chi_sqrt_value, 1)
                self._var_qual = round(-10 * np.log10(chi_prob), 2) \
                    if chi_prob else 10000

        return


###############################################################################
class BaseVarSingleProcess(object):
    """
    simple class to repesent a single BaseVar process.
    """

    # class variable for all instances
    samples_id = None
    total_samples = None
    total_subsamcol = None

    def __init__(self, mpileup_files, out_vcf_file, out_cvg_file,
                regions, options, cmm=None):
        """
        Store input file, options and output file name.

        Parameters
        """
        self.mpileup_files = mpileup_files
        self.out_vcf_file = out_vcf_file
        self.out_cvg_file = out_cvg_file
        self.regions = {}
        self.options = options
        self.cmm = cmm

        # store the region into a dict
        for chrid, start, end in regions:

            if chrid not in self.regions:
                self.regions[chrid] = []

            self.regions[chrid].append([start, end])

        # Cache a batch of mpileup file handle which index by tabix
        self.tb_files = [pysam.TabixFile(f) for f in self.mpileup_files]

        # assignment sample id to the class variable
        if not BaseVarSingleProcess.samples_id:
            # must be called in __init__ function!

            (BaseVarSingleProcess.samples_id,
             BaseVarSingleProcess.total_samples,
             BaseVarSingleProcess.total_subsamcol) = self._load_sample_name()
            sys.stderr.write('[INFO] Finish loading sample name.\n')

    def _load_sample_name(self):
        """
        """
        # load all the samples, 2D array
        sample_id = []
        with open(self.options.samplelistfile) as I:

            samplefiles = []
            for r in I:
                if r.startswith('#'): continue
                samplefiles.append(r.strip().split()[0])

            for f in samplefiles:
                with open(f) as I:
                    sample_id.append([s.strip().split()[0] for s in I])

        total_sample = []
        for s in sample_id:
            total_sample.extend(s)

        # loading subsample if provide
        total_subsamcol = []
        if self.options.subsample:
            subsample = []
            with open(self.options.subsample) as I:
                for r in I:
                    if r.startswith('#'): continue
                    subsample.append(r.strip().split()[0])

            subsample = set(subsample)
            for i, s in enumerate(total_sample):
                if s in subsample:
                    total_subsamcol.append(i)  # get index in total_sample

            total_subsamcol = set(total_subsamcol)

        return sample_id, total_sample, total_subsamcol

    def _close_tabix(self):
        for tb in self.tb_files:
            tb.close()

    def run(self):
        """
        Run the process of calling variant and output
        """
        vcf_header = utils.vcf_header_define()
        with open(self.out_vcf_file, 'w') as VCF, open(self.out_cvg_file, 'w') as CVG:

            CVG.write('\t'.join(['#CHROM', 'POS', 'REF', 'Depth'] + self.cmm.BASE +
                                ['Indel', 'FS', 'Strand_cvg']) + '\n')

            # set header
            VCF.write('\n'.join(vcf_header) + '\n')
            VCF.write('\t'.join(['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t'
                                 'INFO\tFORMAT'] +
                                BaseVarSingleProcess.total_samples) + '\n')

            for chrid, regions in sorted(self.regions.items(), key=lambda x: x[0]):
                # ``regions`` is a 2-D array : [[start1,end1], [start2, end2], ...]
                # fetch the position data from each mpileup files
                # ``iter_tokes`` is a list of iterator for each sample's mpileup
                tmp_region = []
                for p in regions: tmp_region.extend(p)  # covert to 1d-array
                tmp_region = sorted(tmp_region)

                start, end = tmp_region[0], tmp_region[-1]
                iter_tokes = []
                sample_info = []

                for i, tb in enumerate(self.tb_files):
                    try:
                        iter_tokes.append(tb.fetch(chrid, start-1, end))
                    except ValueError:
                        if self.cmm.debug:
                            print >> sys.stderr, ("# [WARMING] Empty region",
                                                  chrid, start-1, end,
                                                  self.mpileup_files[i])
                        iter_tokes.append('')

                # Set iteration marker: 1->iterate; 0->donot
                # iterate or hit the end
                go_iter = [1] * len(iter_tokes)
                for start, end in regions:
                    for position in xrange(start, end + 1):

                        sample_info = [mpileup.fetch_next(iter_tokes[i])
                                       if g else sample_info[i]
                                       for i, g in enumerate(go_iter)]

                        ref_base, sample_base, sample_base_qual, strands, indels = (
                            mpileup.fetch_base_by_position(
                                position,
                                BaseVarSingleProcess.samples_id,
                                sample_info,
                                go_iter,
                                iter_tokes,
                                is_scan_indel=True)
                        )

                        # ignore positions if coverage=0 or ref base is 'N' base
                        if not ref_base or ref_base in ['N', 'n']:
                            continue

                        if BaseVarSingleProcess.total_subsamcol:
                            for k, b in enumerate(sample_base):
                                if k not in BaseVarSingleProcess.total_subsamcol:
                                    # set un-selected bases to be 'N' which
                                    # will be filted
                                    sample_base[k] = 'N'

                        self._out_cvg_file(chrid, position, ref_base, sample_base,
                                           strands, indels, CVG)

                        bt = BaseType(ref_base.upper(), sample_base,
                                      sample_base_qual, cmm=self.cmm)
                        bt.lrt()

                        if len(bt.alt_bases()) > 0:
                            self._out_vcf_line(chrid, position, ref_base,
                                               sample_base, strands, bt, VCF)

        self._close_tabix()

    def _out_cvg_file(self, chrid, position, ref_base, sample_base,
                      strands, indels, out_file_handle):
        # coverage info for each position

        base_depth = {b: 0 for b in self.cmm.BASE}
        for k, b in enumerate(sample_base):

            if self.total_subsamcol and k not in self.total_subsamcol:
                # set un-selected bases to be 'N' which will be filted later
                sample_base[k] = 'N'
                continue

            # ignore all bases('*') which not match ``cmm.BASE``
            if b in base_depth:
                base_depth[b] += 1

        # deal with indels
        indel_dict = {}
        for ind in indels:
            indel_dict[ind] = indel_dict.get(ind, 0) + 1

        indel_string = ','.join(
            [k + ':' + str(v) for k, v in indel_dict.items()]) if indel_dict else '.'

        fs, ref_fwd, ref_rev, alt_fwd, alt_rev = 0, 0, 0, 0, 0
        if sample_base:
            base_sorted = sorted(base_depth.items(),
                                 key=lambda x, y: cmp(x[1], y[1]),
                                 reverse=True)

            b1, b2 = base_sorted[0][0], base_sorted[1][0]
            fs, ref_fwd, ref_rev, alt_fwd, alt_rev = strand_bias(
                ref_base,
                [b1 if b1 != ref_base.upper() else b2],
                sample_base,
                strands
            )

        out_file_handle.write('\t'.join(
            [chrid, str(position), ref_base, str(sum(base_depth.values()))] +
            [str(base_depth[b]) for b in self.cmm.BASE] + [indel_string]) +
                  '\t' + str(fs) + '\t' +
                  ','.join(map(str, [ref_fwd, ref_rev, alt_fwd, alt_rev])) + '\n')

        return

    def _out_vcf_line(self, chrid, position, ref_base, sample_base,
                      strands, bt, out_file_handle):
        #
        alt_gt = {b:'./'+str(k+1) for k,b in enumerate(bt.alt_bases())}
        samples = []

        for k, b in enumerate(sample_base):

            # For sample FORMAT
            if b != 'N':
                # For the base which not in bt.alt_bases()
                if b not in alt_gt: alt_gt[b] = './.'
                gt = '0/.' if b==ref_base.upper() else alt_gt[b]

                samples.append(gt+':'+b+':'+strands[k]+':'+
                               str(round(bt.qual_pvalue[k], 6)))
            else:
                samples.append('./.') ## 'N' base

        # Strand bias by fisher exact test
        # Normally you remove any SNP with FS > 60.0 and an indel with FS > 200.0
        fs, ref_fwd, ref_rev, alt_fwd, alt_rev = strand_bias(
            ref_base, bt.alt_bases(), sample_base, strands)

        # base=>[AF, allele depth]
        af = {b:['%f' % round(bt.depth[b]/float(bt.total_depth), 6),
                 bt.depth[b]] for b in bt.alt_bases()}

        info = {'CM_DP': str(int(bt.total_depth)),
                'CM_AC': ','.join(map(str, [af[b][1] for b in bt.alt_bases()])),
                'CM_AF': ','.join(map(str, [af[b][0] for b in bt.alt_bases()])),
                'CM_EAF': ','.join(map(str, [bt.eaf[b] for b in bt.alt_bases()])),
                'FS': str(fs),
                'SB_REF': str(ref_fwd)+','+str(ref_rev),
                'SB_ALT': str(alt_fwd)+','+str(alt_rev)}

        out_file_handle.write('\t'.join([chrid, str(position), '.', ref_base,
                         ','.join(bt.alt_bases()), str(bt.var_qual()),
                         '.' if bt.var_qual() > self.cmm.QUAL_THRESHOLD else 'LowQual',
                         ';'.join([k+'='+v for k, v in sorted(
                            info.items(), key=lambda x:x[0])]),
                            'GT:AB:SO:BP'] + samples) + '\n')
        return


###############################################################################
class BaseVarMultiProcess(multiprocessing.Process):
    """
    simple class to represent a single BaseVar process, which is run as part of
    a multi-process job.
    """
    def __init__(self, mpileup_files, out_vcf_file, out_cvg_file,
                regions, options, cmm=None):

        """
        Constructor.
        """
        multiprocessing.Process.__init__(self)
        self.single_process = BaseVarSingleProcess(mpileup_files,
                                                   out_vcf_file,
                                                   out_cvg_file,
                                                   regions,
                                                   options,
                                                   cmm=cmm)

    def run(self):
        """ Run the BaseVar process"""
        self.single_process.run()


