# cython: profile=True
"""The codes and functions for defining a GenomeRegion type

Author: shujia Huang
Date: 2020-02-24 10:36:55
"""
from libc.stdio cimport fprintf, stderr
from libc.stdlib cimport exit, EXIT_FAILURE
from libc.stdlib cimport malloc, realloc, free

cdef bint genome_region_array_init(GenomeRegionArray *dest, const size_t size):
    """Allocate a StringArray of size"""
    dest.array = <GenomeRegion *> malloc(size * sizeof(GenomeRegion))
    if dest.array == NULL:
        fprintf(stderr, "[ERROR] Could not allocate memory for GenomeRegionArray")
        exit(EXIT_FAILURE)

    dest.size = 0  # We don't put anything in here yet, just allocate memory
    dest.__capacity = size

    # success
    return True

cdef bint genome_region_init_by_pylist(GenomeRegionArray *dest, list py_array):
    """Init the StringArray by a python_list.
    
    ``py_array``: list like. The element in py_array is [chrom, start, end]
    """
    cdef size_t size = len(py_array)
    genome_region_array_init(dest, size)

    cdef unsigned i = 0
    cdef GenomeRegion temp
    cdef bytes _py_chrom
    for i in range(size):
        _py_chrom = py_array[i][0]
        temp.chrom = _py_chrom
        temp.start = py_array[i][1]
        temp.end = py_array[i][2]
        genome_region_array_append(dest, temp)

    return True

cdef void genome_region_array_destroy(GenomeRegionArray *dest):
    """deallocate the array of String"""
    cdef unsigned int index = 0
    if dest.array != NULL and dest.__capacity > 0:
        free(dest.array)
        dest.size = 0
        dest.__capacity = 0

    return

cdef void genome_region_cpy(GenomeRegion *dest, const GenomeRegion src):
    """Copy the GenomeRegion from 'src' to 'dest' """
    dest.chrom = src.chrom
    dest.start = src.start
    dest.end = src.end
    return

cdef void genome_region_array_append(GenomeRegionArray *dest, const GenomeRegion src):
    """Append a new GenomeRegion into dest array, re-allocate the array if necessary.
    """
    cdef GenomeRegion * temp = NULL
    if dest.size == dest.__capacity:
        temp = <GenomeRegion *> realloc(dest.array, 2 * dest.__capacity * sizeof(GenomeRegion))
        if temp == NULL:
            fprintf(stderr, "[ERROR]Could not re-allocate GenomeRegionArray.")
            exit(EXIT_FAILURE)
        else:
            dest.array = temp
            dest.__capacity *= 2

    genome_region_cpy(&dest.array[dest.size], src)
    dest.size += 1

    return