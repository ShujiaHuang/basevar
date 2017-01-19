"""
"""
from variantcaller import BaseType

if __name__ == '__main__':

    ind_base = ['A', 'A', 'G', 'T', 'G', 'A']
    ind_qual_phredscale = [30, 30, 30, 5, 10, 30]
    bt = BaseType('C', ind_base, ind_qual_phredscale)
    bt.lrt()

    bt.debug()

