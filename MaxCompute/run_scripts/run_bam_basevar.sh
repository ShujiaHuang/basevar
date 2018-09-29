#!/bin/bash

set -e

start_idx=$1
end_idx=$2

INFO_FILE="DONE"
LOG_DIR="./log"
PARTITION_POS_LENTH=6000000

is_continue=$1
if [ "$is_continue" == "--new" ]
then
    rm -f $INFO_FILE
fi

#OSS_INPUT_PATH="genomedata/testdata/140K_BAM/bam"
#OSS_INPUT_PATH="genomedata/1M_CMDB/20180920_PopGen"
OSS_INPUT_PATH="genomedata/1M_CMDB/20180927_PopGen"
#OSS_INPUT_PATH="genomedata/testdata2"
OSS_HOST_INTERNAL="oss-cn-shenzhen-internal.aliyuncs.com"
OSS_HOST="oss-cn-shenzhen.aliyuncs.com"
OSS_ID="LTAIk3YBbHCA8EWk"
OSS_KEY="SGqUx92FF5rVebMDOc3OaZKlWmL811"

#DATASET_NAME="test_14w_time"
#DATASET_NAME="20180920_PopGen"
DATASET_NAME="20180927_PopGen"
FASTA_FILE="Homo_sapiens_assembly38.fasta.gz"
ODPS_CMD="odpscmd --config=/apsarapangu/disk2/tianli.tl/huada/base_var/odps_conf/odps_config.ini.sz"
#chrid="chr8"

$ODPS_CMD -e "
--fs -mkv basevar_ref;
--fs -put /apsarapangu/disk2/tianli.tl/huada/base_var/BaseVar/MaxCompute/run_scripts/log/hg19.fasta.gz basevar_ref/hg19;
--add VOLUMEFILE /basevar_ref/hg19/hg19.fasta.gz as hg19.fasta.gz;
--add VOLUMEFILE /basevar_ref/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta.gz as Homo_sapiens_assembly38.fasta.gz -f;
--add py /apsarapangu/disk2/tianli.tl/huada/base_var/BaseVar/MaxCompute/basevar_part_bam.py -f;
--add jar /tmp/oss-input-1.0.0.jar -f;
--add file /tmp/target_sample.txt -f;
--add file /tmp/popgroup.txt -f;
--add py /apsarapangu/disk2/tianli.tl/huada/base_var/BaseVar/MaxCompute/expand_target_part_bam.py -f;
--create function expand_target_part_bam_0 as expand_target_part_bam.Expand_0 using expand_target_part_bam.py,target_sample.txt -f;
--create function expand_target_part_bam_1 as expand_target_part_bam.Expand_1 using expand_target_part_bam.py,target_sample.txt -f;
--create function expand_target_part_bam_2 as expand_target_part_bam.Expand_2 using expand_target_part_bam.py,target_sample.txt -f;
--add py /apsarapangu/disk2/tianli.tl/huada/base_var/BaseVar/MaxCompute/basetype.py -f;
--add py /apsarapangu/disk2/tianli.tl/huada/base_var/BaseVar/MaxCompute/algorithm.py -f;
--create function basevar_part_bam as basevar_part_bam.BaseVar using basevar_part_bam.py,basetype.py,algorithm.py,scipy.zip,target_sample.txt,popgroup.txt -f;
"

read_input_oss() {
    input_name=$1
    fasta_file=$2
    oss_host=$3
    oss_id=$4
    oss_key=$5
    oss_path=$6
    subdir=$7

    external_table="oss_${input_name}_${subdir}"
    inner_table="${input_name}_partitioned"

    $ODPS_CMD -e "
-- Defind external table
DROP TABLE IF EXISTS $external_table;
CREATE EXTERNAL TABLE IF NOT EXISTS $external_table
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
    'fasta_file' = '${fasta_file}'
    --'sample_filter' = 'target_sample.txt' #ignore sample_filter
)
LOCATION 'oss://${oss_id}:${oss_key}@${oss_host}/${oss_path}/${subdir}'
USING 'oss-input-1.0.0.jar,${fasta_file}';
"

    $ODPS_CMD -e "
CREATE TABLE IF NOT EXISTS $inner_table
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
    chr STRING,
    shard STRING
)
"

    $ODPS_CMD -e "
set odps.sql.mapper.memory=6144;
--set odps.sql.mapper.memory=10240;
--set odps.sql.udf.jvm.memory=10240;
set odps.sql.udf.jvm.memory=8192;
set odps.sql.reshuffle.dynamicpt=false;
set odps.sql.planner.mode=lot;
--set odps.sql.unstructured.data.split.size=32;
INSERT INTO TABLE ${inner_table} PARTITION (chr, shard)
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
    chrid chr,
    CAST(pos / $PARTITION_POS_LENTH AS BIGINT) shard
FROM
    ${external_table}
WHERE
    base_ref <> 'N';

DROP TABLE IF EXISTS ${external_table};
"
}


expand_pos_data() {
    input_name=$1
    chr_id=$2
    shard_id=$3
    expand_name="${input_name}_expand"
    inner_table="${input_name}_partitioned"
    #expand
    $ODPS_CMD -e "
set odps.sql.map.aggr=false;
--set odps.sql.udf.timeout=3600;
set odps.sql.executionengine.batch.rowcount=16;
set odps.sql.planner.mode=lot;
--set odps.sql.runtime.mode=ganjiang;
set odps.sql.mapper.split.size=1280;
--set odps.sql.reducer.instances=6;
set odps.sql.udf.python.memory=4096;
set odps.sql.reducer.memory=6144;
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
    chr STRING,
    shard STRING
);

INSERT OVERWRITE TABLE ${expand_name} PARTITION (chr='$chr_id', shard='$shard_id')
SELECT
    chrid,
    pos,
    base_ref,
    expand_target_part_bam_0(sample_name, read_base, read_quality, mapping_quality, read_pos_rank, indel, strand) AS part0,
    expand_target_part_bam_1(sample_name, read_base, read_quality, mapping_quality, read_pos_rank, indel, strand) AS part1,
    expand_target_part_bam_2(sample_name, read_base, read_quality, mapping_quality, read_pos_rank, indel, strand) AS part2
FROM
    $inner_table
WHERE
    chr = '$chr_id'
AND
    shard = '$shard_id'
GROUP BY
    chrid, pos, base_ref;
"
}

do_both() {
    input_name=$1
    chr_id=$2
    shard_id=$3
    expand_name="${input_name}_expand"
    res_name="${input_name}_result"
    $ODPS_CMD -e "
create table if not exists $res_name (line string) PARTITIONED BY (chr STRING, shard STRING, method STRING);
set odps.sql.planner.mode=lot;
insert overwrite table $vcf_name PARTITION (chr='$chr_id', shard='$shard_id', method)
select
    basevar_part_bam('both', chrid, pos, base_ref, part0, part1, part2) as (line, method)
from
    $expand_name
where
    chr = '$chr_id'
AND
    shard = '$shard_id'
"
}

do_vcf() {
    input_name=$1
    chr_id=$2
    shard_id=$3
    vcf_name="${input_name}_vcf"
    expand_name="${input_name}_expand"
    $ODPS_CMD -e "
create table if not exists $vcf_name (line string) PARTITIONED BY (chr STRING, shard STRING);
set odps.sql.planner.mode=lot;
insert overwrite table $vcf_name PARTITION (chr='$chr_id', shard='$shard_id')
select
    basevar_part_bam('vcf', chrid, pos, base_ref, part0, part1, part2) as vcf
from
    $expand_name
where
    chr = '$chr_id'
AND
    shard = '$shard_id'
"
}

do_cvg() {
    input_name=$1
    chr_id=$2
    cvg_name="${input_name}_cvg"
    expand_name="${input_name}_expand"
    $ODPS_CMD -e "
create table if not exists $cvg_name (line string) PARTITIONED BY (chr STRING);
set odps.sql.planner.mode=lot;
insert overwrite table $cvg_name PARTITION (chr='$chr_id')
select
    basevar_part_bam('coverage', chrid, pos, base_ref, part0, part1, part2) as cvg
from
    $expand_name
where
    chr = '$chr_id'
"
}

for((i=$start_idx;i<=$end_idx;i++))
do
    subdir="sample_list_$i"
    echo "$start_idx,$end_idx,read_input_oss:$subdir:$(date)" >> time_record
    grep -e "^$subdir\$" $INFO_FILE || (read_input_oss ${DATASET_NAME} $FASTA_FILE $OSS_HOST_INTERNAL $OSS_ID $OSS_KEY ${OSS_INPUT_PATH} $subdir && echo $subdir >> $INFO_FILE)
done

exit

for subdir in $(osscmd --host=${OSS_HOST} --id=${OSS_ID} --key=${OSS_KEY} listalldir oss://${OSS_INPUT_PATH} 2>/dev/null | grep -v ": ")
do
    echo "read_input_oss:$subdir:$(date)" >> time_record
    grep $subdir $INFO_FILE || (read_input_oss ${DATASET_NAME} $FASTA_FILE $OSS_HOST_INTERNAL $OSS_ID $OSS_KEY ${OSS_INPUT_PATH} $subdir && echo $subdir >> $INFO_FILE)
done

for chr in $($ODPS_CMD -e "show partitions ${DATASET_NAME}_partitioned" 2>/dev/null | grep chr | awk '{print substr($0, 5)}')
do
    echo "expand_pos_data:$chr:$(date)" >> time_record
    expand_pos_data $DATASET_NAME $chr
done

for chr in $($ODPS_CMD -e "show partitions ${DATASET_NAME}_partitioned" 2>/dev/null | grep chr | awk '{print substr($0, 5)}')
do
    echo "do_vcf:$chr:$(date)" >> time_record
    do_vcf $DATASET_NAME $chr
done

for chr in $($ODPS_CMD -e "show partitions ${DATASET_NAME}_partitioned" 2>/dev/null | grep chr | awk '{print substr($0, 5)}')
do
    echo "do_cvg:$chr:$(date)" >> time_record
    do_cvg $DATASET_NAME $chr
done

for chr in $($ODPS_CMD -e "show partitions ${DATASET_NAME}_partitioned" 2>/dev/null | grep chr | awk '{print substr($0, 5)}')
do
    echo "basevar_exttable_vcf:$chr:$(date)" >> time_record
    ./basevar_exttable.sh ${DATASET_NAME}_vcf $chr
done

for chr in $($ODPS_CMD -e "show partitions ${DATASET_NAME}_partitioned" 2>/dev/null | grep chr | awk '{print substr($0, 5)}')
do
    echo "basevar_exttable_cvg:$chr:$(date)" >> time_record
    ./basevar_exttable.sh ${DATASET_NAME}_cvg $chr
done

echo "alldone:$(date)" >> time_record
exit


# External table
$ODPS_CMD -e "
DROP TABLE IF EXISTS oss_${DATASET_NAME}_input;
CREATE EXTERNAL TABLE IF NOT EXISTS oss_${DATASET_NAME}_input
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
    'fasta_file' = '$FASTA_FILE',
    'sample_filter' = 'target_sample.txt'
)
LOCATION 'oss://LTAIk3YBbHCA8EWk:SGqUx92FF5rVebMDOc3OaZKlWmL811@oss-cn-shenzhen-internal.aliyuncs.com/genomedata/testdata/fusion_test/'
-- LOCATION 'oss://LTAIzpgTEbfsEote:QAnDU4tSlMF1ewsvjLh7w04WFCsFE0@oss-cn-shenzhen-internal.aliyuncs.com/genomedata-sj/nifty_bwa/'
--LOCATION 'oss://LTAIk3YBbHCA8EWk:SGqUx92FF5rVebMDOc3OaZKlWmL811@oss-cn-shenzhen-internal.aliyuncs.com/nifty-140k/bamfile/'
USING 'oss-input-1.0.0.jar,$FASTA_FILE,target_sample.txt';
"

# # Inner table && partition by chrid
$ODPS_CMD -e "
CREATE TABLE IF NOT EXISTS ${DATASET_NAME}_partitioned
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
"
#
#

$ODPS_CMD -e "
set odps.sql.udf.timeout=1200;
set odps.sql.mapper.memory=6144;
set odps.sql.udf.jvm.memory=10240;
set odps.sql.reshuffle.dynamicpt=false;
set odps.sql.planner.mode=lot;
set odps.sql.unstructured.data.split.size=3072;
INSERT OVERWRITE TABLE ${DATASET_NAME}_partitioned PARTITION (chr)
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
    oss_${DATASET_NAME}_input
WHERE
    base_ref <> 'N'
"

 #expand
$ODPS_CMD -e "
set odps.sql.map.aggr=false;
--set odps.sql.udf.timeout=3600;
set odps.sql.reducer.instances=6;
set odps.sql.executionengine.batch.rowcount=16;
set odps.sql.planner.mode=lot;
--set odps.sql.runtime.mode=ganjiang;
set odps.sql.mapper.split.size=1280;
set odps.sql.udf.python.memory=4096;
set odps.sql.reducer.memory=6144;
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

--INSERT OVERWRITE TABLE ${DATASET_NAME}_expand PARTITION (chr='$chrid')
INSERT OVERWRITE TABLE ${DATASET_NAME}_expand PARTITION (chr='0')
SELECT
    chrid,
    pos,
    base_ref,
    expand_target_part_bam_0(sample_name, read_base, read_quality, mapping_quality, read_pos_rank, indel, strand) AS part0,
    expand_target_part_bam_1(sample_name, read_base, read_quality, mapping_quality, read_pos_rank, indel, strand) AS part1,
    expand_target_part_bam_2(sample_name, read_base, read_quality, mapping_quality, read_pos_rank, indel, strand) AS part2
FROM
    ${DATASET_NAME}_partitioned
--WHERE
--    chr = '$chrid'
GROUP BY
    chrid, pos, base_ref;
"
# #
# $ODPS_CMD -e "
# create table if not exists ${DATASET_NAME}_cvg (line string) PARTITIONED BY (chr STRING);
# set odps.sql.planner.mode=lot;
# insert overwrite table ${DATASET_NAME}_cvg PARTITION (chr='$chrid')
# select
#     basevar_part_bam('coverage', chrid, pos, base_ref, part0, part1, part2) as cvg
# from
#     ${DATASET_NAME}_expand
# where
#     chr = '$chrid'
# "
# ./basevar_exttable.sh ${DATASET_NAME}_cvg $chrid
#
$ODPS_CMD -e "
create table if not exists ${DATASET_NAME}_vcf (line string) PARTITIONED BY (chr STRING);
set odps.sql.planner.mode=lot;
--insert overwrite table ${DATASET_NAME}_vcf PARTITION (chr='$chrid')
insert overwrite table ${DATASET_NAME}_vcf PARTITION (chr='0')
select
    basevar_part_bam('vcf', chrid, pos, base_ref, part0, part1, part2) as vcf
from
    ${DATASET_NAME}_expand
--where
--    chr = '$chrid'
"
#./basevar_exttable.sh ${DATASET_NAME}_vcf $chrid
./basevar_exttable.sh ${DATASET_NAME}_vcf 0

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
