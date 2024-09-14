# cython: profile=True
"""The codes and functions for defining a GenomeRegion type

Author: shujia Huang
Date: 2020-02-24 10:36:55
"""
from libc.stdio cimport fprintf, sprintf, stderr
from libc.stdlib cimport exit, EXIT_FAILURE
from libc.stdlib cimport calloc, realloc, free

cdef char *make_region_str(const GenomeRegion src):
    """return the GenomeRegion to be a string: chrom:start-end"""
    cdef unsigned int MAX_SIZE = 200
    cdef char *region = <char *>calloc(MAX_SIZE, sizeof(char))
    if region == NULL:
        fprintf(stderr, "[ERROR] Memory allocate failure in `make_region_str(const GenomeRegion src)`\n")
        exit(EXIT_FAILURE)

    sprintf(region, "%s:%ld-%ld", src.chrom, src.start, src.end)

    # trailing space
    # cdef unsigned int size = MAX_SIZE
    # while(size > 0 and isspace(region[size-1])):
    #     size -=1
    # region[size] = '\0'
    return region

cdef bint genome_region_array_init(GenomeRegionArray *dest, const size_t size):
    """Allocate a StringArray of size"""
    dest.array = <GenomeRegion *> calloc(size, sizeof(GenomeRegion))
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