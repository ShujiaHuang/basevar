"""header for utils.pyx"""

cdef void fast_merge_files(list temp_file_names, basestring final_file_name, bint is_del_raw_file)
cdef list generate_regions_by_process_num(list regions, int process_num, bint convert_to_2d)