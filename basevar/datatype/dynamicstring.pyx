"""Code for dynamic string.

Author: Shujia Huang
Date: 2020-02-28 15:48:44

"""
from libc.stdio cimport fprintf, sprintf, stderr
from libc.stdlib cimport exit, EXIT_FAILURE
from libc.stdlib cimport calloc, realloc, free
from libc.string cimport memcpy, strlen

cdef int _ds_resize(dstring *dest, size_t size):
    cdef char *temp
    if dest.__capacity < size+1:
        temp = <char *>realloc(dest.s, 2 * (size+1) * sizeof(char))
        if temp == NULL:
            fprintf(stderr, "[ERROR]Could not re-allocate dstring.\n")
            exit(EXIT_FAILURE)

        dest.s = temp
        dest.__capacity = 2 * (size + 1)

    return 0

cdef int dstring_init(dstring *dest, const size_t size):
    """Allocate a dstring by size"""
    dest.s = <char *>calloc(size+1, sizeof(char))  # +1 for the '\0' in C string
    if dest.s == NULL:
        fprintf(stderr, "[ERROR] Could not allocate memory when dstring_init.\n")
        exit(EXIT_FAILURE)

    dest.size = 0  # We don't put anything in here yet, just allocate memory
    dest.__capacity = size
    return 0

cdef int dstring_destroy(dstring *dest):
    dest.size = 0
    dest.__capacity = 0
    free(dest.s)
    return 0

cdef int dstring_strcpy(dstring *dest, const char *src):
    cdef size_t size = strlen(src)
    if dest.__capacity > 0 or dest.s != NULL:
        dstring_destroy(dest)

    dstring_init(dest, size)
    dstring_append(dest, src)
    return 0

cdef int dstring_append(dstring *dest, const char *src):
    """Append c string (char*)
    """
    cdef size_t src_size = strlen(src)
    if dest.__capacity < dest.size + src_size + 1:
        _ds_resize(dest, dest.size + src_size + 1)

    memcpy(dest.s + dest.size, src, src_size)
    dest.size += src_size
    dest.s[dest.size] = '\0'  # Mark the end of the c string

    return 0

cdef int dstring_append_char(dstring *dest, char c):
    if dest.__capacity < dest.size + 1 + 1:
        _ds_resize(dest, dest.size + 1 + 1)

    dest.s[dest.size] = c
    dest.size += 1
    dest.s[dest.size] = '\0'
    return 0

cdef int dstring_append_int(dstring *dest, int c):
    cdef char buf[16]
    sprintf(buf, "%d", c)
    dstring_append(dest, buf)
    return 0

cdef int dstring_append_long(dstring *dest, long int c):
    cdef char buf[32]
    sprintf(buf, "%ld", c)
    dstring_append(dest, buf)
    return 0
