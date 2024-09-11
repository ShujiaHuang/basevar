#!/bin/bash
python create_pipeline.py -R human_reference.fa --ref_fai human_reference.fa.fai -c chr1 --delta 5000000 -t 20 -L input_bamfile.list -o output_directory/ > run_chr1_basevar_shell.sh

# Add population groups
python create_pipeline.py -R human_reference.fa --ref_fai human_reference.fa.fai -c chr1 --delta 5000000 -t 20 -L input_bamfile.list -G pop_group.list -o output_chr1_directory/ > run_basevar.sh


