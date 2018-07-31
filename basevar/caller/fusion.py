"""
This is a module for creating fusion file by BAM/CRAM
"""

import sys
import pysam


class FusionElement(object):
    def __init__(self, chrid, mapq, strand_orientation):
        self.chrid = chrid
        self.mapq = mapq
        self.strand_orientation = strand_orientation
        self.start = 0
        self.end = 0
        self.alt = ''  # <NON_REF>, [.], [ACGT]
        self.read_first_position = 0  # record the first_position
        self.base_quality = ''


class Fusion(object):
    """
    Create fusion file for each alignement file
    """

    def __init__(self, in_ref_file, in_align_file):
        """
        Init setting.

        Parameters
        ==========

            ``in_ref_file``: String.
                Reference fasta format.

            ``in_align_file``: String
                Alignement files, BAM/CRAM format

        """
        self.ref_file_hd = pysam.FastaFile(in_ref_file)
        self.align_file_hd = None
        # if in_align_file.endswith('.bam'):
        #     self.align_file_hd = pysam.AlignmentFile(in_align_file, 'rb')
        #
        # elif in_align_file.endswith('.cram'):
        #     self.align_file_hd = pysam.AlignmentFile(in_align_file, 'rc')

        if in_align_file.endswith('.bam') or in_align_file.endswith('.cram'):
            self.align_file_hd = pysam.AlignmentFile(in_align_file)
        else:
            sys.stderr.write('[ERROR] Input file: %s is not BAM nor CRAM.\n' % in_align_file)
            self._close_file()
            sys.exit(1)

    def generate_fusion(self):

        start_postion = {}
        for read in self.align_file_hd:

            # Ignore the read which mapping quality lower than 30
            if read.mapq < 30:
                continue

            start_postion[read.reference_name] = start_postion.get(read.reference_name, 0)

            # get the reference sequence
            fa = self.ref_file_hd.fetch(read.reference_name)
            for fusion in self.get_fusion_region(fa, read, start_postion[read.reference_name]):
                yield fusion

            if read.reference_end > start_postion[read.reference_name]:
                start_postion[read.reference_name] = read.reference_end

        self._close_file()

    def get_fusion_region(self, fa, read, start_postion):

        """
        Get fusion region by scan the alignment reads

        Parameter
        =========

            ``fa``: String
                Fasta sequence

            ``read``:  pysam.libcalignedsegment.AlignedSegment
                Aligned reads

            ``start_position``: Integer
                The position which start to find aligment info from ``read``

        Return
        ======
            return an array content fusion-element

        """
        fusion_region = []
        fusion = FusionElement(read.reference_name, read.mapq,
                               '+' if not read.is_reverse else '-')

        """
        The cigar string order in the array is "MIDNSHP=X" followed by a
        field for the NM tag. If the NM tag is not present, this
        field will always be 0.

            +-----+--------------+-----+
            |M    |BAM_CMATCH    |0    |
            +-----+--------------+-----+
            |I    |BAM_CINS      |1    |
            +-----+--------------+-----+
            |D    |BAM_CDEL      |2    |
            +-----+--------------+-----+
            |N    |BAM_CREF_SKIP |3    |
            +-----+--------------+-----+
            |S    |BAM_CSOFT_CLIP|4    |
            +-----+--------------+-----+
            |H    |BAM_CHARD_CLIP|5    |
            +-----+--------------+-----+
            |P    |BAM_CPAD      |6    |
            +-----+--------------+-----+
            |=    |BAM_CEQUAL    |7    |
            +-----+--------------+-----+
            |X    |BAM_CDIFF     |8    |
            +-----+--------------+-----+
            |B    |BAM_CBACK     |9    |
            +-----+--------------+-----+
            |NM   |NM tag        |10   |
            +-----+--------------+-----+

        If no cigar string is present, empty arrays will be archived.
        """
        clip_head, clip_tail = 0, 0

        # CLIP can only happen in read head or read tail.
        if read.cigar[0][0] == 4: # 4:SOFT_CLIP Head
            clip_head = read.cigar[0][1]

        if len(read.cigar) > 1 and read.cigar[-1][0] == 4: # 4:SOFT_CLIE Tail
            clip_tail = read.cigar[-1][1]

        is_first = True
        is_indel = False

        ref_pos_not_in_insertion = 0
        # loop the aligne-pairs to find the alignment status of read and reference.
        # qpos and ref_pos are all 0-base system
        for qpos, ref_pos in read.aligned_pairs[clip_head:len(read.aligned_pairs)-clip_tail]:

            if ref_pos is not None:
                ref_pos_not_in_insertion = ref_pos

            # Ignore the position which has been scanned
            if ref_pos_not_in_insertion < start_postion:
                continue

            ref_base = fa[ref_pos].upper() if ref_pos is not None else '.'
            read_base = read.seq[qpos].upper() if qpos is not None else '.'
            base_qual = read.qual[qpos] if qpos is not None else '!'  # '!'-33 is 0 for indel

            # Initial Indel fusion element
            if (qpos is None or ref_pos is None) and not is_indel:
                fusion = FusionElement(read.reference_name, read.mapq,
                                       '+' if not read.is_reverse else '-')

            if qpos is None:  # deletion

                is_indel = True

                # 1-base system: fusion.start = fusion_region[-1].end + 1
                fusion.start = fusion_region[-1].end  # set to be 0-base
                fusion.end = ref_pos + 1
                fusion.alt = '.'
                fusion.base_quality = base_qual

            elif ref_pos is None:  # insertion

                is_indel = True

                # stay 0-base coordinate system
                fusion.start = fusion_region[-1].end - 1
                fusion.end = fusion_region[-1].end
                fusion.alt += read_base
                fusion.base_quality += base_qual

            else:  ## Matching

                if is_indel:

                    # store Indel information
                    fusion_region.append(fusion)
                    fusion = FusionElement(read.reference_name, read.mapq,
                                           '+' if not read.is_reverse else '-')

                    # Refresh if it's the indel region in last one of fusion_region
                    is_first = True

                is_indel = False

                # data setting
                fusion.start = ref_pos  # stay 0-base coordinate system
                fusion.end = ref_pos + 1  # End position 0-base => 1-base system

                # just record the first rank position in read, we can easily
                # get other rank when the position is continuous
                fusion.read_first_position = qpos + 1
                fusion.base_quality = base_qual

                if read_base != ref_base:  # candidate variants

                    fusion.alt = read_base
                else:
                    # Although read_base == ref_base, it may still be
                    # reference base or not.
                    fusion.alt = '<NON_REF>'

            if is_indel:
                continue

            if not is_first:

                if (fusion_region[-1].alt == '<NON_REF>') and (read_base == ref_base) and \
                        (fusion_region[-1].base_quality == base_qual):
                    # extending the NON_REF Region
                    fusion_region[-1].end = fusion.end

                else:

                    fusion_region.append(fusion)
                    fusion = FusionElement(read.reference_name, read.mapq,
                                           '+' if not read.is_reverse else '-')

            else:
                is_first = False

                fusion_region.append(fusion)
                fusion = FusionElement(read.reference_name, read.mapq,
                                       '+' if not read.is_reverse else '-')

        return fusion_region

    def _close_file(self):
        self.ref_file_hd.close()
        self.align_file_hd.close()


if __name__ == '__main__':

    # Just for testing the Fusion module
    callfusion = Fusion(sys.argv[1], sys.argv[2])
    for fusion in callfusion.generate_fusion():
        print '\t'.join(map(str, [fusion.chrid,
                                  fusion.start,
                                  fusion.end,
                                  fusion.alt,
                                  fusion.mapq,
                                  fusion.strand_orientation,
                                  fusion.read_first_position,
                                  fusion.base_quality]
                            )
                        )
