## 1 资源文件准备
这里首先把要用到的资源文件、函数添加到odps project中，如果后续没有改动，这一步不用重复执行（每次数据的target_sample.list文件是不一样的，需要更新）。

```
fs -mkv basevar_ref;
fs -put /tmp/Homo_sapiens_assembly38.fasta.gz basevar_ref/Homo_sapiens_assembly38;
add VOLUMEFILE /basevar_ref/Homo_sapiens_assembly38 as Homo_sapiens_assembly38.fasta.gz;
add py /tmp/basevar_part_bam.py -f;
add jar /tmp/oss-input-1.0.0.jar -f;
add file /tmp/target_sample.list -f;
add py /tmp/expand_target_part_bam.py;
create function expand_target_part_bam as expand_target_part_bam.Expand using expand_target_part_bam.py,target_sample.list;
add py /tmp/basetype.py -f;
add py /tmp/algorithm.py -f;
create function basevar_part_bam as basevar_part_bam.BaseVar using basevar_part_bam.py,basetype.py,algorithm.py,scipy.zip;
```

## 2 通过外表读取oss数据，并进行预处理pileup

第一步从oss中读取数据并存入odps内表。由于数据量大，读取时按chrid拆分，方便后续作业按染色体分别运行。读取过程中同时进行预处理，完成pileup。

```
-- 创建外表，指定列名、基因序列参考文件(fasta_file)、sample文件(target_sample.list)和oss数据所在路径信息
CREATE EXTERNAL TABLE IF NOT EXISTS oss_nifty_input
(
    sample_name STRING,
    chrid STRING,
    pos STRING,
    base_ref STRING,
    read_base STRING,
    read_quality STRING,
    mapping_quality STRING,
    read_pos_rank STRING,
    indel STRING,
    strand STRING
)
STORED BY 'com.aliyun.odps.exttable.bam.BamHandler'
WITH SERDEPROPERTIES
(
    'fasta_file' = 'Homo_sapiens_assembly38.fasta.gz',
    'sample_filter' = 'target_sample.list'
)
LOCATION 'oss://LTAIzpgTEbfsEote:QAnDU4tSlMF1ewsvjLh7w04WFCsFE0@oss-cn-shenzhen-internal.aliyuncs.com/genomedata-sj/nifty_bwa/'
USING 'oss-input-1.0.0.jar,Homo_sapiens_assembly38.fasta.gz,target_sample.list';

-- 创建内部分区表，与外表对应
CREATE TABLE IF NOT EXISTS nifty_partitioned
(
    sample_name STRING,
    chrid STRING,
    pos STRING,
    base_ref STRING,
    read_base STRING,
    read_quality STRING,
    mapping_quality STRING,
    read_pos_rank STRING,
    indel STRING,
    strand STRING
)
PARTITIONED BY
(
    chr STRING
)

-- 外表导入数据到分区表，按chrid分区
set odps.sql.udf.timeout=1200;
set odps.sql.mapper.memory=6144;
set odps.sql.udf.jvm.memory=10240;
set odps.sql.reshuffle.dynamicpt=false;
set odps.sql.planner.mode=lot;
set odps.sql.unstructured.data.split.size=3072;
INSERT OVERWRITE TABLE nifty_partitioned PARTITION (chr)
SELECT
    sample_name,
    chrid,
    pos,
    base_ref,
    read_base,
    read_quality,
    mapping_quality,
    read_pos_rank,
    indel,
    strand,
    chrid chr
FROM
    oss_nifty_input
WHERE
    base_ref <> 'N'
```

## 3 Expand 成宽表

上一步数据读取进来后，要按照位点把所有的sample数据汇集到一个字段里，用于计算vcf和cvg。

这里我们直接按照位点做group by，然后用UDF把sample的数据拼接合并起来。由于拼接后数据量较大，odps对单个字段的大小有限制，所以这里分成3列来存储。

```
-- 建表，宽表数据量大，一列存不下，拆分成3列 part0 / part1 / part2
set odps.sql.map.aggr=false;
set odps.sql.reducer.instances=3000;
set odps.sql.executionengine.batch.rowcount=16;
set odps.sql.planner.mode=lot;
set odps.sql.mapper.split.size=1280;
set odps.sql.udf.python.memory=4096;
set odps.sql.reducer.memory=6144;
CREATE TABLE IF NOT EXISTS nifty_expand
(
    chrid STRING,
    pos STRING,
    base_ref STRING,
    part0 STRING,
    part1 STRING,
    part2 STRING
)
PARTITIONED BY
(
    chr STRING
);

-- 逐条染色体合并到宽表 (指定chr为 chr1 / chr2 / ... / chr22 / chrX)
INSERT OVERWRITE TABLE nifty_expand PARTITION (chr='chr1')
SELECT
    chrid,
    pos,
    base_ref,
    expand_target_part_bam(sample_name, read_base, read_quality, mapping_quality, read_pos_rank, indel, strand, 0) AS part0,
    expand_target_part_bam(sample_name, read_base, read_quality, mapping_quality, read_pos_rank, indel, strand, 1) AS part1,
    expand_target_part_bam(sample_name, read_base, read_quality, mapping_quality, read_pos_rank, indel, strand, 2) AS part2
FROM
    nifty_partitioned
WHERE
    chr = 'chr1'
GROUP BY
    chrid, pos, base_ref;
```

## 4 计算 vcf 和 cvg
由于上一步我们把数据拆成了3列，这里计算之前就要先把3列合并起来，然后用华大同学提供的算法代码来计算最终的结果。所以基于原先的basevar，调整成basevar_part_bam，进行最终结果计算输出。

```
-- 建vcf结果表
create table if not exists nifty_vcf (line string) PARTITIONED BY (chr STRING);
set odps.sql.planner.mode=lot;
-- 逐条染色体计算vcf (指定chr为 chr1 / chr2 / ... / chr22 / chrX)
insert overwrite table nifty_vcf PARTITION (chr='chr1')
select
    basevar_part_bam('vcf', chrid, pos, base_ref, part0, part1, part2) as vcf
from
    nifty_expand
where
    chr = 'chr1'
```

```
-- 建cvg结果表
create table if not exists nifty_cvg (line string) PARTITIONED BY (chr STRING);
set odps.sql.planner.mode=lot;
-- 逐条染色体计算cvg (指定chr为 chr1 / chr2 / ... / chr22 / chrX)
insert overwrite table nifty_cvg PARTITION (chr='chr1')
select
    basevar_part_bam('coverage', chrid, pos, base_ref, part0, part1, part2) as cvg
from
    nifty_expand
where
    chr = 'chr1'
```

