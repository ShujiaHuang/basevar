"""Open general files

Author: Shujia Huang
Date: 2019-06-08 02:13:26

"""
import os
import gzip
import heapq

from basevar.io.BGZF.bgzf import BGZFile

def _expanded_open(path, mode):
    try:
        return open(path, mode)
    except IOError:
        return open(os.path.expanduser(path), mode)


def Open(file_name, mode, compress_level=9, isbgz=True):
    """
    Function that allows transparent usage of dictzip, gzip and
    ordinary files
    """
    if file_name.endswith(".gz") or file_name.endswith(".GZ"):
        file_dir = os.path.dirname(file_name)
        if not os.path.exists(file_dir):
            file_name = os.path.expanduser(file_name)

        return BGZFile(file_name, mode) if isbgz else gzip.GzipFile(file_name, mode, compress_level)
    else:
        return _expanded_open(file_name, mode)


class FileForQueueing(object):
    def __init__(self, the_file, line, is_del_raw_file=False):
        """
        Store the file, and init current value
        """
        self.the_file = the_file
        self.finishedReadingFile = False
        self.is_del_raw_file = is_del_raw_file
        self.heap = []

        line = line
        cols = line.strip().split()
        chrom = cols[0]

        # Where possible, convert chromosome names into
        # integers for sorting. If not possible, use
        # original names.
        try:
            chrom = int(chrom.upper().strip("CHR"))
        except Exception:
            pass

        pos = int(cols[1])
        heapq.heappush(self.heap, (chrom, pos, line))

        while not self.finishedReadingFile and len(self.heap) < 100:

            try:
                line = self.the_file.next()
                cols = line.strip().split()
                chrom = cols[0]

                try:
                    chrom = int(chrom.upper().strip("CHR"))
                except Exception:
                    pass

                pos = int(cols[1])
            except StopIteration:
                self.finishedReadingFile = True
                break

            heapq.heappush(self.heap, (chrom, pos, line))

        # take the top line
        self.chrom, self.pos, self.line = heapq.heappop(self.heap)

    def __cmp__(self, other):
        """
        Comparison function. Utilises the comparison function defined in
        the AlignedRead class.
        """
        def _comparision(a, b):
            if a < b:
                return -1
            elif a > b:
                return 1
            else:
                return 0

        return _comparision(self.chrom, other.chrom) or _comparision(self.pos, other.pos)

    def __del__(self):
        """
        Destructor
        """
        self.the_file.close()

        if self.is_del_raw_file:
            os.remove(self.the_file.name)

    def next(self):
        """
        Increment the iterator and yield the new value. Also, store the
        current value for use in the comparison function.
        """
        if not self.finishedReadingFile:

            try:
                line = self.the_file.next()
                cols = line.strip().split()
                chrom = cols[0]

                # Where possible, convert chromosome names into
                # integers for sorting. If not possible, use
                # original names.
                try:
                    chrom = int(chrom.upper().strip("CHR"))
                except Exception:
                    pass

                pos = int(cols[1])
                heapq.heappush(self.heap, (chrom, pos, line))

            except StopIteration:
                self.finishedReadingFile = True

        if len(self.heap) != 0:
            # Now take the top line
            self.chrom, self.pos, self.line = heapq.heappop(self.heap)
        else:
            raise StopIteration