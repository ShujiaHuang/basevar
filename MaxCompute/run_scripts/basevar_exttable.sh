#!/bin/bash

set -e

table_name=$1  #basevar_140k_sample_chr11_cvg
chr_name=$2
shard_id=$3
method=$4

real_table=bgi_max_sz.$table_name

ODPS_CMD='odpscmd --config=/apsarapangu/disk2/tianli.tl/huada/base_var/odps_conf/odps_config.ini.sz'

#osscmd --host=oss-cn-shenzhen.aliyuncs.com --id=LTAIk3YBbHCA8EWk --key=SGqUx92FF5rVebMDOc3OaZKlWmL811 mkdir oss://nifty-140k/basevar_results/${table_name}_${chr_name}
#osscmd --host=oss-cn-shenzhen.aliyuncs.com --id=LTAIzpgTEbfsEote --key=QAnDU4tSlMF1ewsvjLh7w04WFCsFE0 mkdir oss://genomedata-sj/test_data/fusiontest_results/${table_name}_${chr_name}
osscmd --host=oss-cn-shenzhen.aliyuncs.com --id=LTAIk3YBbHCA8EWk --key=SGqUx92FF5rVebMDOc3OaZKlWmL811 mkdir oss://genomedata/test_data/fusiontest_results/${table_name}/${chr_name}

pos_table="${table_name}_${chr_name}_${shard_id}_${method}_withpos"
$ODPS_CMD -e "
DROP TABLE IF EXISTS ${pos_table};
CREATE TABLE IF NOT EXISTS ${pos_table} AS
SELECT 
    CAST(split_part(line, '\t', 2, 2) AS BIGINT) pos, 
    line 
FROM 
    ${real_table}
WHERE
    chr='$chr_name'
AND
    shard='$shard_id'
AND
    method='$method'
"

min_max=$($ODPS_CMD -e "select min(pos), max(pos) from ${pos_table}" | grep '[0-9]' | grep -v '_')
min_pos=$(echo $min_max | awk '{print $2}')
max_pos=$(echo $min_max | awk '{print $4}')
echo $min_max
echo $min_pos $max_pos
table_size=$($ODPS_CMD -e "desc $real_table" | grep Size | awk '{print $6}')
echo $table_size
per_partition=$((2 * 1024 * 1024 * 1024))
partition_count=$(($table_size / $per_partition + 1))
interval=$((($max_pos - $min_pos) / $partition_count + 1))

echo $partition_count
echo $interval
#exit

partitioned_pos_table="${pos_table}_partitioned"
$ODPS_CMD -e "
DROP TABLE IF EXISTS ${partitioned_pos_table};
CREATE TABLE IF NOT EXISTS ${partitioned_pos_table} AS
SELECT 
    CAST((pos - $min_pos) / $interval AS BIGINT) partition_id,
    pos,
    line
FROM
    ${pos_table};
DROP TABLE IF EXISTS ${pos_table};
"

ossout_table="${partitioned_pos_table}_ossout"
$ODPS_CMD -e "
DROP TABLE IF EXISTS ${ossout_table};
CREATE EXTERNAL TABLE IF NOT EXISTS ${ossout_table} (
    partition_id BIGINT,
    idx          BIGINT,
    line         STRING
)
STORED BY 'com.aliyun.odps.exttable.handler.BasevarHandler'
--LOCATION 'oss://LTAIk3YBbHCA8EWk:SGqUx92FF5rVebMDOc3OaZKlWmL811@oss-cn-shenzhen-internal.aliyuncs.com/nifty-140k/basevar_results/${table_name}_${chr_name}/'
--LOCATION 'oss://LTAIzpgTEbfsEote:QAnDU4tSlMF1ewsvjLh7w04WFCsFE0@oss-cn-shenzhen-internal.aliyuncs.com/genomedata-sj/basevar_results/${table_name}_${chr_name}/'
LOCATION 'oss://LTAIk3YBbHCA8EWk:SGqUx92FF5rVebMDOc3OaZKlWmL811@oss-cn-shenzhen-internal.aliyuncs.com/genomedata/test_data/fusiontest_results/${table_name}/${chr_name}/'
USING 'oss-input-1.0.0.jar';
"

$ODPS_CMD -e "
set odps.sql.executionengine.batch.rowcount=16;
set odps.sql.reducer.instances=$partition_count;
insert overwrite table ${ossout_table}
select * from ${partitioned_pos_table}
distribute by partition_id
sort by partition_id, pos;
"

echo "minmax:" $min_pos $max_pos ",table_size:" $table_size ",partition_count:" $partition_count ",interval:" $interval

