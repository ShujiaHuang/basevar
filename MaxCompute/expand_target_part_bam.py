#coding:utf-8
from odps.udf import annotate
from odps.udf import BaseUDAF
#from odps.distcache import get_cache_file
import sys

@annotate('string,string,string,string,string,string,string->string')
class Expand(BaseUDAF):
    def __init__(self, part_idx):
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
        self.part_idx = part_idx
        self.buffer_size = max(0, min(self.sample_count - self.part_size * self.part_idx, self.part_size))
        #print >> sys.stderr, 'Expand init size: %s, part_idx: %s' % (self.buffer_size, self.part_idx)

    def new_buffer(self):
        #print >> sys.stderr, 'new_buffer size: %s, part_idx: %s' % (self.buffer_size, self.part_idx)
        return [''] * self.buffer_size


    def iterate(self, buffer, sample_name, read_base, read_quality, mapping_quality, read_pos_rank, indel, strand):
        #print >> sys.stderr, 'iterate'
        idx = self.name_idx.get(sample_name, -1)
        idx_of_part = idx - self.part_size * self.part_idx
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
        if self.buffer_size <= 0:
            return ''
        for i, s in enumerate(buffer):
            if not s:
                buffer[i] = 'N\t0\t0\t0\t\t.' # default value
        return '\t'.join(buffer[0:self.buffer_size])


@annotate('string,string,string,string,string,string,string->string')
class Expand_0(Expand):
    def __init__(self):
        Expand.__init__(self, 0)

@annotate('string,string,string,string,string,string,string->string')
class Expand_1(Expand):
    def __init__(self):
        Expand.__init__(self, 1)


@annotate('string,string,string,string,string,string,string->string')
class Expand_2(Expand):
    def __init__(self):
        Expand.__init__(self, 2)
