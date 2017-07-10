# MaxCompute BaseVar README

## Prepare OSS

upload all mpileup gz files and sample name files into one oss folder,
eg. huada-test/nifty/
generate a split.txt file and upload it into the same oss folder.

## Create external table in MaxCompute using OSS files

merge all sample name files into one: `all_names.txt`
and add it as a maxcompute resource file:
`add file all_names.txt;`

compile mpileip_extractor: `mvn package`
and upload mpileup-extractor-1.0.jar as a maxcompute resource jar:
`add jar mpileup-extractor-1.0.jar;`

create external table:
```
set odps.task.major.version=2dot0_demo_flighting;
set odps.sql.planner.mode=lot;
set odps.sql.ddl.odps2=true;
set odps.sql.preparse.odps2=lot;

CREATE EXTERNAL TABLE IF NOT EXISTS oss_nifty
(
chrid string,
pos string,
base_ref string,
sample_name string,
c1 string,
c2 string,
c3 string
)
STORED BY 'com.aliyun.odps.poc.MpileupHandler'
WITH SERDEPROPERTIES
(
'odps.properties.rolearn'='acs:ram::1817850323806830:role/aliyunodpsdefaultrole',
'oss.user.defined.files.splits'='split.txt'
)
LOCATION 'oss://oss-cn-beijing-internal.aliyuncs.com/huada-test/nifty/'
USING 'mpileup-extractor-1.0.jar';
```

# for debug

```
CREATE EXTERNAL TABLE IF NOT EXISTS oss_small
(
chrid string,
pos string,
base_ref string,
sample_name string,
c1 string,
c2 string,
c3 string
)
STORED BY 'com.aliyun.odps.poc.MpileupHandler'
WITH SERDEPROPERTIES
(
'odps.properties.rolearn'='acs:ram::1817850323806830:role/aliyunodpsdefaultrole',
'oss.user.defined.files.splits'='split_small.txt'
)
LOCATION 'oss://oss-cn-beijing-internal.aliyuncs.com/huada-test/nifty/'
USING 'mpileup-extractor-1.0.jar';
```
