"""
"""
cdef int call_vcfconcat(list args):
    cdef int num = len(args)
    cdef char **c_args = NULL
    c_args =  <char**>(calloc(num, sizeof(char*)))

    args = map(bytes, args)
    cdef int i = 0, size = 0
    for i in range(num):
        size = len(args[i])
        c_args[i] = <char*>(calloc(size, sizeof(char)))
        c_args[i] = args[i]

    cdef int flag = main_vcfconcat(num, c_args)

    if c_args != NULL:
        for i in range(num):
            if c_args[i] != NULL:
                free(c_args[i])
                c_args[i] = NULL

        free(c_args)
        c_args = NULL

    return flag

def main(args):
    return call_vcfconcat(args)


