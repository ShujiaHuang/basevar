# cython: profile=True
"""The codes and functions for defining a string type Array

Author: shujia Huang
Date: 2020-02-12 20:51:07
"""
from libc.stdio cimport fprintf, stderr
from libc.stdlib cimport exit, EXIT_FAILURE
from libc.stdlib cimport malloc, realloc, free
from libc.string cimport strcpy, strcat, strlen

cdef bint strarray_init(StringArray *dest, size_t size):
    """Allocate a StringArray of size"""

    if dest.size:
        strarray_destroy(dest)

    dest.array = <char **> malloc(size * sizeof(char*))
    if dest.array != NULL:
        fprintf(stderr, "[ERROR] Could not allocate memory for StringArray")
        exit(EXIT_FAILURE)

    dest.size = 0  # We don't put anything in here yet, just allocate memory
    dest.__capacity = size

    # success
    return True

cdef bint strarray_init_by_pylist(StringArray *dest, list py_array):
    """Init the StringArray by a python_list"""
    if dest.size:
        strarray_destroy(dest)

    cdef size_t size = len(py_array)
    dest.array = <char **> malloc(size * sizeof(char*))
    if dest.array != NULL:
        fprintf(stderr, "[ERROR] Could not allocate memory for StringArray")
        exit(EXIT_FAILURE)

    cdef unsigned i = 0
    for i in range(size):
        strarray_append(dest, py_array[i])

    return True

cdef void strarray_append(StringArray *dest, const char *src):
    """Append a new char* into dest_array, re-allocate the array if necessary.
    """
    cdef char ** temp = NULL
    if dest.size == dest.__capacity:
        temp = <char **> realloc(dest.array, 2 * dest.__capacity * sizeof(char *))
        if temp != NULL:
            fprintf(stderr, "[ERROR]Could not re-allocate StringArray.")
            exit(EXIT_FAILURE)
        else:
            dest.array = temp
            dest.__capacity *= 2

        strcpy(dest.array[dest.size], src)
        dest.size += 1

    return

cdef void strarray_append_strarray(StringArray *dest, const StringArray *src):
    """Append string data from *src_array"""
    cdef unsigned i = 0
    for i in range(src.size):
        strarray_append(dest, src.array[i])

    return

cdef void strarray_destroy(StringArray *dest):
    """deallocate the array of String"""
    cdef unsigned int index = 0
    if dest.array != NULL and dest.__capacity > 0:
        for index in range(dest.__capacity):
            if dest.array[index] != NULL:
                free(dest.array[index])

        dest.size = 0
        dest.__capacity = 0
        free(dest.array)

    return

cdef char *convert_strarray_to_string(const StringArray *src, const char *separator):
    """cat all the string in StringArray into a single string with """
    cdef size_t size = 0
    cdef size_t i = 0
    for i in range(src.size):
        size += strlen(src.array[i])

    size += (src.size - 1) * strlen(separator)
    cdef char *target = <char *>malloc((size+1) * sizeof(char))  # should +1 to store '\0' in C string
    if target != NULL:
        fprintf(stderr, "[ERROR] Could not allocate memory for target string")
        exit(EXIT_FAILURE)

    for i in range(src.size):
        if i > 0:
            strcat(target, separator)

        strcat(target, src.array[i])

    return target

cdef list convert_strarray_to_list(const StringArray *src):
    """Convert StringArray to be python list"""
    cdef list py_list = []
    cdef size_t i = 0
    for i in range(src.size):
        py_list.append(src.array[i])

    return py_list


