#!/bin/bash

set -e

DATASET_NAME="nifty_140k_bam"
FASTA_FILE="hg19.fasta.gz"
ODPS_CMD="odpscmd --config=/apsarapangu/disk2/tianli.tl/huada/base_var/odps_conf/odps_config.ini.sz"
chrid="chr1"

# $ODPS_CMD -e "
# fs -mkv basevar_ref;
# fs -put /apsarapangu/disk2/tianli.tl/huada/base_var/BaseVar/MaxCompute/run_scripts/log/hg19.fasta.gz basevar_ref/hg19;
# add VOLUMEFILE /basevar_ref/hg19/hg19.fasta.gz as hg19.fasta.gz;
# add py /apsarapangu/disk2/tianli.tl/huada/base_var/BaseVar/MaxCompute/run_scripts/odps_resource/basevar_part_bam.py -f;
# add jar /tmp/oss-input-1.0.0.jar -f;
# add file /tmp/target_sample.list -f;
# add py /apsarapangu/disk2/tianli.tl/huada/base_var/BaseVar/MaxCompute/run_scripts/odps_resource/expand_target_part_bam.py;
# create function expand_target_part_bam as expand_target_part_bam.Expand using expand_target_part_bam.py,target_sample.list;
# add py /apsarapangu/disk2/tianli.tl/huada/base_var/BaseVar/MaxCompute/basetype.py -f;
# add py /apsarapangu/disk2/tianli.tl/huada/base_var/BaseVar/MaxCompute/algorithm.py -f;
# create function basevar_part_bam as basevar_part_bam.BaseVar using basevar_part_bam.py,basetype.py,algorithm.py,scipy.zip;
# "

# External table
# $ODPS_CMD -e "
# DROP TABLE IF EXISTS oss_${DATASET_NAME}_input;
# CREATE EXTERNAL TABLE IF NOT EXISTS oss_${DATASET_NAME}_input
# (
#     sample_name STRING,
#     chrid STRING,
#     pos STRING,
#     base_ref STRING,
#     read_base STRING,
#     read_quality STRING,
#     mapping_quality STRING,
#     read_pos_rank STRING,
#     indel STRING,
#     strand STRING
# )
# STORED BY 'com.aliyun.odps.exttable.bam.BamHandler'
# WITH SERDEPROPERTIES
# (
#     'fasta_file' = '$FASTA_FILE'
# )
# LOCATION 'oss://LTAIk3YBbHCA8EWk:SGqUx92FF5rVebMDOc3OaZKlWmL811@oss-cn-shenzhen-internal.aliyuncs.com/nifty-140k/bamfile/'
# USING 'oss-input-1.0.0.jar,$FASTA_FILE';
# "
#
# # Inner table && partition by chrid
# $ODPS_CMD -e "
# CREATE TABLE IF NOT EXISTS ${DATASET_NAME}_partitioned
# (
#     sample_name STRING,
#     chrid STRING,
#     pos STRING,
#     base_ref STRING,
#     read_base STRING,
#     read_quality STRING,
#     mapping_quality STRING,
#     read_pos_rank STRING,
#     indel STRING,
#     strand STRING
# )
# PARTITIONED BY
# (
#     chr STRING
# )
# "
#
# $ODPS_CMD -e "
# set odps.sql.udf.timeout=1200;
# set odps.sql.mapper.memory=6144;
# set odps.sql.udf.jvm.memory=10240;
# set odps.sql.reshuffle.dynamicpt=false;
# set odps.sql.planner.mode=lot;
# set odps.sql.unstructured.data.split.size=2048;
# INSERT OVERWRITE TABLE ${DATASET_NAME}_partitioned PARTITION (chr)
# SELECT
#     sample_name,
#     chrid,
#     pos,
#     base_ref,
#     read_base,
#     read_quality,
#     mapping_quality,
#     read_pos_rank,
#     indel,
#     strand,
#     chrid chr
# FROM
#     oss_${DATASET_NAME}_input
# WHERE
#     base_ref <> 'N'
# "

# expand
$ODPS_CMD -e "
set odps.sql.map.aggr=false;
set odps.sql.udf.timeout=3600;
set odps.sql.reducer.instances=3000;
set odps.sql.executionengine.batch.rowcount=1;
CREATE TABLE IF NOT EXISTS ${DATASET_NAME}_expand
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

INSERT OVERWRITE TABLE ${DATASET_NAME}_expand PARTITION (chr='$chrid')
SELECT
    chrid,
    pos,
    base_ref,
    expand_target_part_bam(sample_name, read_base, read_quality, mapping_quality, read_pos_rank, indel, strand, 0) AS part0,
    expand_target_part_bam(sample_name, read_base, read_quality, mapping_quality, read_pos_rank, indel, strand, 1) AS part1,
    expand_target_part_bam(sample_name, read_base, read_quality, mapping_quality, read_pos_rank, indel, strand, 2) AS part2
FROM
    ${DATASET_NAME}_partitioned
WHERE
    chr = '$chrid'
GROUP BY
    chrid, pos, base_ref;
"

$ODPS_CMD -e "create table if not exists ${DATASET_NAME}_cvg (line string) PARTITIONED BY (chr STRING);
set odps.sql.planner.mode=lot;
insert overwrite table ${DATASET_NAME}_cvg PARTITION (chr='$chrid')
select
    basevar_part_bam('coverage', chrid, pos, base_ref, part0, part1, part2) as cvg
from
    ${DATASET_NAME}_expand
where
    chr = '$chrid'
"
./basevar_exttable.sh ${DATASET_NAME}_cvg $chrid

$ODPS_CMD -e "create table if not exists ${DATASET_NAME}_vcf (line string) PARTITIONED BY (chr STRING);
set odps.sql.planner.mode=lot;
insert overwrite table ${DATASET_NAME}_vcf PARTITION (chr='$chrid')
select
    basevar_part_bam('vcf', chrid, pos, base_ref, part0, part1, part2) as vcf
from
    ${DATASET_NAME}_expand
where
    chr = '$chrid'
"
./basevar_exttable.sh ${DATASET_NAME}_vcf $chrid

#$ODPS_CMD -e "
#CREATE EXTERNAL TABLE IF NOT EXISTS oss_nifty_140k
#(
#chrid string,
#pos string,
#base_ref string,
#sample_name string,
#c1 string,
#c2 string,
#c3 string
#)
#STORED BY 'com.aliyun.odps.exttable.mpileup.MpileupHandler'
#LOCATION 'oss://LTAIzpgTEbfsEote:QAnDU4tSlMF1ewsvjLh7w04WFCsFE0@oss-cn-shenzhen-internal.aliyuncs.com/genomedata-sj/nifty_140k/mpileup/'
#USING 'oss-input-1.0.0.jar';"
