#coding:utf-8
from odps.udf import annotate
from odps.udf import BaseUDAF
from odps.distcache import get_cache_file

@annotate('string,string,string,string->string')
class Expand(BaseUDAF):

    def __init__(self):
        res_file = get_cache_file('all_names.txt')
        name_idx = dict()
        i = 0
        for line in res_file:
            name_idx[line.strip()] = i
            i += 1
        res_file.close()
        self.name_idx = name_idx
        self.total_name = i + 1

    def new_buffer(self):
        return [""] * self.total_name

    def iterate(self, buffer, sample_name, c1, c2, c3):
        s = "\t".join([c1, c2, c3])
        buffer[self.name_idx[sample_name]] = s

    def merge(self, buffer, pbuffer):
        for i, s in enumerate(pbuffer):
            if s:
                if buffer[i]:
                    raise Exception('merge buffer error: %s' % i)
                buffer[i] = s

    def terminate(self, buffer):
        for i, s in enumerate(buffer):
            if not s:
                buffer[i] = "N\t!\t." # default value
        return "\t".join(buffer)
