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

## Expand OSS external table to a MaxCompute wide table

add python resource file and create function:
```
add py expand_udf.py;
create function expand as expand_udf.Expand using expand_udf.py,all_names.txt;
```

create wide table:
```
create table nifty_expand as select chrid, pos, base_ref, expand(sample_name, c1, c2, c3) as one from nifty group by chrid, pos, base_ref;
```

## BaseVar Calculation

prepare scipy environment for MaxCompute, download scipy whl package first
`wget http://mirrors.aliyun.com/pypi/packages/ae/94/28ca6f9311e2351bb68da41ff8c1bc8f82bb82791f2ecd34efa953e60576/scipy-0.19.0-cp27-cp27m-manylinux1_x86_64.whl#md5=0e49f7fc8d31c1c79f0a4d63b29e8a1f`

add scipy whl package as a archive resource:
`add archive scipy-0.19.1-cp27-cp27m-manylinux1_x86_64.whl as scipy.zip;`

create udf to for calculating BaseVar coverage:
```
add py basevar.py;
add py basetype.py;
add py mpileup.py;
add py algorithm.py;
create function basevar as basevar.BaseVar using mpileup.py,basevar.py,basetype.py,algorithm.py,scipy.zip;
```

run query:
```
set odps.task.major.version=2dot0_demo_flighting;
set odps.sql.planner.mode=lot;
set odps.sql.ddl.odps2=true;
set odps.sql.preparse.odps2=lot;

set odps.pypy.enabled=false;
set odps.isolation.session.enable = true;

create table nifty_cvg (line string);
create table nifty_vcf (line string);

from nifty_expand
insert overwrite table nifty_cvg
select basevar('coverage', chrid, pos, base_ref, one) as line
insert overwrite table nifty_vcf
select basevar('vcf', chrid, pos, base_ref, one) as line;
```

or generate coverage/vcf table separately:
```
create table nifty_cvg as select basevar('coverage', chrid, pos, base_ref, one) as cvg from nifty_expand;

create table nifty_vcf as select basevar('vcf', chrid, pos, base_ref, one) as vcf from nifty_expand;
```
