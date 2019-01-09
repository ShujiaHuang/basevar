
## 从OSS下载参考序列
```bash
oss://genomedata/human_reference/hg19_NC_012920/hg19.NC_012920.fasta
oss://genomedata/human_reference/hg19_NC_012920/hg19.NC_012920.fasta.fai
```

## 从BAM生成Fusion，并使用bgzip压缩，tabix建索引

```bash
OIFS=$IFS
for i in `cat bam100.list` ; do

    IFS='/' read -r -a ADD <<< "$i"
    IFS='.' read -r -a sample <<< "${ADD[13]}"

    sample_id=${sample[0]}
    sub_out_dir="/Volumes/Macintosh_HD/Users/huangshujia/PycharmProjects/BaseVar/tests/data/140k_thalassemia_brca_bam/fusion/${ADD[11]}/${ADD[12]}"

done
```


## Basetype

```bash
time python ../../../basevar/BaseVar.py basetype -I bam100/00alzqq6jw.bam -I bam100/09t3r9n2rg.bam -I bam100/0fkpl1p55b.bam -I bam100/13dg1gvsfk.bam -I bam100/17phildszl.bam -I bam100/1dbpgqt0dq.bam -I bam100/1kyws27hoc.bam -I bam100/1ych8rmufr.bam -I bam100/4e56w6ezsx.bam -I bam100/51rwla2fps.bam -R ../hg19.NC_012920.fasta --regions chr11:5246595-5248428,chr17:41197764-41276135 --batch-count 10 --output-vcf test.vcf --output-cvg test.cvg.tsv --pop-group sample_group.info --nCPU 4 -L bam90.list
```

## Just create batch file

```bash
time python ../../../basevar/BaseVar.py basetype -I bam100/00alzqq6jw.bam -I bam100/09t3r9n2rg.bam -I bam100/0fkpl1p55b.bam -I bam100/13dg1gvsfk.bam -I bam100/17phildszl.bam -I bam100/1dbpgqt0dq.bam -I bam100/1kyws27hoc.bam -I bam100/1ych8rmufr.bam -I bam100/4e56w6ezsx.bam -I bam100/51rwla2fps.bam -R ../hg19.NC_012920.fasta --regions chr11:5246595-5248428,chr17:41197764-41276135 --batch-count 10 --output-batch-file test.batchfile.gz --nCPU 4 -L bam90.list
```

## Basetype from batchfile
```bash
time python ../../../basevar/BaseVar.py basetypebatch -I test.batchfile.gz -R ../hg19.NC_012920.fasta --regions chr11:5246595-5248428,chr17:41197764-41276135 --output-vcf test.vcf --output-cvg test.cvg.tsv --pop-group sample_group.info --nCPU 4
```
