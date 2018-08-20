#coding:utf-8
from odps.udf import annotate
from odps.udf import BaseUDAF
#from odps.distcache import get_cache_file

@annotate('string,string,string,string,string,string,string,bigint->string')
class Expand(BaseUDAF):

    def __init__(self):
        from odps.distcache import get_cache_file
        res_file = get_cache_file('target_sample.txt')
        self.name_idx = dict()
        i = 0
        for line in res_file:
            self.name_idx[line.strip()] = i
            i += 1
        res_file.close()
        self.sample_count = i
        self.part_size = 50000
        self.part_idx = -1
        self.buffer_size = min(self.sample_count, self.part_size)

    def new_buffer(self):
        return [''] * self.buffer_size

    def iterate(self, buffer, sample_name, read_base, read_quality, mapping_quality, read_pos_rank, indel, strand, part_idx):
        self.part_idx = part_idx
        idx = self.name_idx.get(sample_name, -1)
        idx_of_part = idx - self.part_size * part_idx
        if idx_of_part >= 0 and idx_of_part < self.buffer_size:
            s = '\t'.join([read_base, read_quality, mapping_quality, read_pos_rank, indel, strand])
            buffer[idx_of_part] = s

    def merge(self, buffer, pbuffer):
        for i, s in enumerate(pbuffer):
            if s:
                if buffer[i]:
                    raise Exception('merge buffer error: %s' % i)
                buffer[i] = s

    def terminate(self, buffer):
        buffer_length = min(self.sample_count - self.part_size * self.part_idx, self.part_size)
        if buffer_length <= 0:
            return ''
        for i, s in enumerate(buffer):
            if not s:
                buffer[i] = 'N\t0\t0\t0\t\t.' # default value
        return '\t'.join(buffer[0:buffer_length])
