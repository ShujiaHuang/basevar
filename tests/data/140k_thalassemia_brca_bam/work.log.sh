../../../bin/basevar basetype --mapq=10 --min-af=0.05 --batch-count=1 --thread=1 --regions=CHROMOSOME_I:900-1200 --pop-group=sample_group.info --output-vcf vz.vcf --output-cvg t.cvg -R ce.fa.gz -I range.bam -I range.bam > log2

#chr11:5246595-5248428,chr13:32890633-32972781,chr16:222869-227506,chr17:41197764-41276135
../../../bin/basevar basetype --mapq=10 --min-af=0.05 --batch-count=20 --thread=4 --regions=chr11:5246595-5248428,chr17:41197764-41276135 --pop-group=sample_group.info --output-vcf tt.vcf.gz --output-cvg tt.cvg.gz -R ~/Projects/BaseVar/tests/data/hg19.NC_012920.fasta.gz -I bam100/00alzqq6jw.bam -I bam100/09t3r9n2rg.bam -I bam100/0fkpl1p55b.bam -I bam100/13dg1gvsfk.bam -I bam100/17phildszl.bam -I bam100/1dbpgqt0dq.bam -I bam100/1kyws27hoc.bam -I bam100/1ych8rmufr.bam -I bam100/4e56w6ezsx.bam -I bam100/51rwla2fps.bam -L bam90.list > log2



