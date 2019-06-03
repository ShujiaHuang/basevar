"""Header for basetype.pyx
"""
cdef class BaseType:
    cdef:
        bytes _ref_base
        bytes[:] _alt_bases
        float total_depth
        double _var_qual
        double min_af
        double[:,:] ind_allele_likelihood
        list qual_pvalue
        dict af_by_lrt
        dict depth
        void _set_init_ind_allele_likelihood(self, bytes[:] bases, char[:] base_element)

