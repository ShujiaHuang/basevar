"""
A class for output VCF file. PyVCF does not able to add or update information 
fields for sample's FORMAT field. That make us have to create another classes
(like this) to handle that problem
"""
import re

class Header(object):

    def __init__(self, hInfo = None): 
        """
        VCF header information
        """
        self.header = {}
        if hInfo and (type(hInfo) is not dict): 
            raise ValueError ('The data type should be "dict" in class '
                              'of "VCFHeader", but found %s' % str(type(hInfo)))
        if hInfo:
            self.header = hInfo
        
    def add(self, mark, id, num, type, description):
        key = '##%s=<ID=%s' % (mark, id)
        val = ('##%s=<ID=%s,Number=%s,Type=%s,Description="%s">' % 
              (mark, id, num if num is not None else '.', type, description))

        self.header[key] = val
        return self

    def record(self, headline):
        if re.search (r'^##fileformat', headline):
            tag = '###'
        elif re.search (r'^#CHROM', headline):
            tag = '#CHROM'
        else:
            tag = headline.split(',')[0]

        self.header[tag] = headline


class Info(object): 

    def __init__(self, info = None):
        """
        INOF fields information
        """
        self.info = {}
        if info and (type(info) is not dict): 
            raise ValueError ('The data type should be "dict" in class '
                              'of "VCFInfo", but found %s' % str(type(info)))
        if info: self.info = info

    def add(self, key, context):
        self.info[key] = context

        return self


class Context(object): 

    def __init__(self):

        """
        VCF comtext
        """
        self.chrom = None
        self.pos   = None  
        self.Id    = None
        self.ref   = None  
        self.alt   = []
        self.qual  = None    
        self.filter = []
        self.info   = {} 
        self.format = None 
        self.sample = []

    def print_context(self):
        """
        """
        if self.chrom:

            print ('\t'.join([self.chrom,
                     str(self.pos),
                     '.' if not self.Id else self.Id,
                     self.ref,
                     ','.join(self.alt),
                     str(self.qual),
                     '.' if not self.filter else ','.join(self.filter),
                     '.' if not self.info else ';'.join(v
                            for v in sorted(self.info.values())),
                     ':'.join(self.format),
                     '\t'.join(self.sample)]))


def calcuInbreedCoeff(gt):
    """
    Calculating the inbreeding coefficient by GT fields of VCF.

    Args:
        `gt`: A list. Genotype fields of all the samples.
    """
    ref_count, het_count, hom_count, n = 0, 0, 0, 0
    for g in gt:
        gs = g.split('/') if '/' in g else g.split('|')
        if '.' not in g: n += 1
        if '.' in g:
            # Do nothing
            pass

        elif g == '0/0' or g == '0|0':
            # Reference
            ref_count += 1

        elif gs[0] == gs[1]:
            # homo
            hom_count += 1
        else:
            # hete
            het_count += 1

    if n == 0: n = 1
    p = (2.0 * ref_count + het_count) / (2.0 * n) # expected REF allele freq
    q = 1.0 - p # expected alternative allele frequency
    # Inbreeding coefficient: the het_count VS expected of het_count
    expected_het_count = 2.0 * p * q * n if p * q > 0 else 1
    inbf = 1.0 - het_count / expected_het_count
    
    return inbf

def loadPedigree(pedigree_file):
    """
    """
    if not pedigree_file: return {} # 

    pedigree = {}
    for line in open(pedigree_file):
        # 1006 1006-05 1006-01 1006-02 0 0
        col = line.strip('\n').split()
        if col[1] in pedigree:
            raise ValueError('[ERROR] %s is already in "pedigree". Your file ' 
                             'may have the duplication sample name.' % col[1])
        pedigree[col[1]] = [col[2], col[3]]

    return pedigree



