#ifndef B_STREAM_H
#define B_STREAM_H

#include "htslib/kseq.h"
// fastq parsing
// KSEQ_INIT(gzFile, gzread)
KSEQ_INIT(BGZF *, bgzf_read)
//KSTREAM_INIT( gzFile, gzread, 16384)

#endif
