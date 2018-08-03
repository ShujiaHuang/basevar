
## 从BAM生成Fusion，并使用bgzip压缩，tabix建索引

```bash
OIFS=$IFS
for i in `cat bam100.list` ; do

    IFS='/' read -r -a ADD <<< "$i"
    IFS='.' read -r -a sample <<< "${ADD[13]}"

    sample_id=${sample[0]}
    sub_out_dir="/Volumes/Macintosh_HD/Users/huangshujia/PycharmProjects/BaseVar/tests/data/140k_thalassemia_brca_bam/fusion/${ADD[11]}/${ADD[12]}"

    if [ ! -d $sub_out_dir ]
        then mkdir -p $sub_out_dir
    fi

    # echo "$sub_out_dir/${sample_id}.fusion.bed.gz"

    echo "** dealing $sample_id **" && \
    time python ../../../basevar/BaseVar.py fusion -R hg19.NC_012920.fasta -I $i -O $sub_out_dir/${sample_id}.fusion.bed && \
    bgzip -f $sub_out_dir/${sample_id}.fusion.bed && \
    tabix -f -p bed $sub_out_dir/${sample_id}.fusion.bed.gz && \
    echo "** done **" 

done
```

## 从Fusion生成VCF（FusionBaseType 测试）

```bash
time python ../../../basevar/BaseVar.py fusionbasetype --fusion-file-list f100.list --pop-group sample_group.info --reference hg19.NC_012920.fasta --region chr11:5246595-5248428,chr17:41197764-41276135 --outprefix testfusion_group 2> log2
```
