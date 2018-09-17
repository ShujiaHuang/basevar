#!/bin/bash

set -e

table_name=$1  #basevar_140k_sample_chr11_cvg
chr_name=$2

real_table=bgi_max_sz.$table_name

ODPS_CMD='odpscmd --config=/apsarapangu/disk2/tianli.tl/huada/base_var/odps_conf/odps_config.ini.sz'

#osscmd --host=oss-cn-shenzhen.aliyuncs.com --id=LTAIk3YBbHCA8EWk --key=SGqUx92FF5rVebMDOc3OaZKlWmL811 mkdir oss://nifty-140k/basevar_results/${table_name}_${chr_name}
#osscmd --host=oss-cn-shenzhen.aliyuncs.com --id=LTAIzpgTEbfsEote --key=QAnDU4tSlMF1ewsvjLh7w04WFCsFE0 mkdir oss://genomedata-sj/test_data/fusiontest_results/${table_name}_${chr_name}
osscmd --host=oss-cn-shenzhen.aliyuncs.com --id=LTAIk3YBbHCA8EWk --key=SGqUx92FF5rVebMDOc3OaZKlWmL811 mkdir oss://genomedata/test_data/fusiontest_results/${table_name}/${chr_name}

$ODPS_CMD -e "
DROP TABLE IF EXISTS ${table_name}_${chr_name}_withpos;
CREATE TABLE IF NOT EXISTS ${table_name}_${chr_name}_withpos AS
SELECT 
    CAST(split_part(line, '\t', 2, 2) AS BIGINT) pos, 
    line 
FROM 
    ${real_table}
WHERE
    chr='$chr_name'
"

min_max=$($ODPS_CMD -e "select min(pos), max(pos) from ${table_name}_${chr_name}_withpos" | grep '[0-9]' | grep -v '_')
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

$ODPS_CMD -e "
DROP TABLE IF EXISTS ${table_name}_${chr_name}_foross;
CREATE TABLE IF NOT EXISTS ${table_name}_${chr_name}_foross AS
SELECT 
    CAST((pos - $min_pos) / $interval AS BIGINT) partition_id,
    pos,
    line
FROM
    ${table_name}_${chr_name}_withpos
"

$ODPS_CMD -e "
DROP TABLE IF EXISTS ${table_name}_${chr_name}_ossout;
CREATE EXTERNAL TABLE IF NOT EXISTS ${table_name}_${chr_name}_ossout (
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
insert overwrite table ${table_name}_${chr_name}_ossout
select * from ${table_name}_${chr_name}_foross
distribute by partition_id
sort by partition_id, pos;
"

echo "minmax:" $min_pos $max_pos ",table_size:" $table_size ",partition_count:" $partition_count ",interval:" $interval

