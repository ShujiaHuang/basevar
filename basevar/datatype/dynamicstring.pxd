"""Code for dynamic string
Author: Shujia Huang
Date: 2020-02-28 15:48:44
"""
# dstring is a simple dynamic string for BaseVar
ctypedef struct dstring:
    size_t size
    size_t __capacity
    char *s

cdef int dstring_init(dstring *dest, const size_t size)
cdef int dstring_destroy(dstring *dest)
cdef int dstring_strcpy(dstring *dest, const char *src)
cdef int dstring_append(dstring *dest, const char *src)
cdef int dstring_append_char(dstring *dest, char c)
cdef int dstring_append_int(dstring *dest, int c)
cdef int dstring_append_long(dstring *dest, long int c)
