"""The codes and functions for defining a string type Array

Author: shujia Huang
Date: 2020-02-12 20:51:07
"""
ctypedef struct StringArray:
    size_t size
    size_t __capacity
    char ** array

cdef bint strarray_init(StringArray *array, size_t size)
cdef bint strarray_init_by_pylist(StringArray *dest, list py_array)
cdef void strarray_append(StringArray *array, const char *data)
cdef void strarray_append_strarray(StringArray *dest, const StringArray *src)
cdef void strarray_destroy(StringArray *array)
cdef char *convert_strarray_to_string(const StringArray *src, const char *delimiter)
cdef list convert_strarray_to_list(const StringArray *src)
