"""The codes and functions for defining a GenomeRegion and GenomeRegionArray type

Author: shujia Huang
Date: 2020-02-24 10:36:55
"""
ctypedef struct GenomeRegion:
    char *chrom
    unsigned long start
    unsigned long end

ctypedef struct GenomeRegionArray:
    size_t size
    size_t __capacity
    GenomeRegion *array

cdef bint genome_region_array_init(GenomeRegionArray *dest, const size_t size)
cdef bint genome_region_init_by_pylist(GenomeRegionArray *dest, list py_array)
cdef void genome_region_array_destroy(GenomeRegionArray *dest)
cdef void genome_region_cpy(GenomeRegion *dest, const GenomeRegion src)
cdef void genome_region_array_append(GenomeRegionArray *dest, const GenomeRegion src)
