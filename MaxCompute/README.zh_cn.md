# MaxCompute BaseVar 算法实现说明文档

## 协作开发环境

为进行 BaseVar 算法的开发，已在 MaxCompute 开通名为 huada_test 的 project。必要的账号、设置记录于开发机的 /home/work/USER/ruibo/huada_test.ini 文件，相应的 odpscmd 也已配置妥当，可以直接使用。

文中代码均已提交至 BaseVar 代码库的 maxcompute 分支，后不再赘述。

## 数据准备

BaseVar 算法的输入是一个具体基因点位，100 万样本的 base/qual/strand 值。

而 vcf 计算完毕并合并之后，是若干 mpileup 文件及配套的 sample name 文件。

sample name 文件较简单，内容示意如下：

```
nl241syo66
xk54r1w51s
...
```

mpileup 文件和 sample name 总是成对出现，内容示意如下：

| 位置              | sample 1, 对应 name nl241syo66 | sample 2, 对应 name xk54r1w51s | mpileup 文件中通常有 200 个 sample |
| --------------- | ---------------------------- | ---------------------------- | --------------------------- |
| chr11 5246595 C | 0 * *                        | 1 , B                        | ...                         |
| ...             | ...                          | ...                          | ...                         |

为了填补算法输入要求和 OSS 文件格式之间的差异，算法的起点是将 OSS 文件转换成要求的表结构。整个数据准备过程可以大致描述为两个步骤：

1. 通过创建 OSS 外表（ [参考资料](https://help.aliyun.com/document_detail/45389.html)）的方式，将 mpileup 的 200 列表转化为单列表，表结构如下：

| chrid | pos     | base_ref | sample_name | b    | q    | s    |
| ----- | ------- | -------- | ----------- | ---- | ---- | ---- |
| chr11 | 5246595 | C        | nl241syo66  | 0    | *    | *    |
| chr11 | 5246595 | C        | xk54r1w51s  | 1    | ,    | B    |
| ...   |         |          |             |      |      |      |

2. 通过 SQL 语句将上述外表聚合为 100w 列的宽表。

### 创建 OSS 外表

首先，将 mpileup/sample name 文件全部集中在  OSS 的**一个目录**下。根据文件名构造类似如下的 split.txt 文件（内容中文件名冒号后面的数字为 magic number，沿用例子里的数字即可），并上传至相同的 OSS 目录下。

```
140k_sample.1.list.samplename:0
140k_sample_test.1.mpileup.gz:67108864
140k_sample.2.list.samplename:0
140k_sample_test.2.mpileup.gz:67108864
140k_sample.3.list.samplename:0
140k_sample_test.3.mpileup.gz:67108864
140k_sample.4.list.samplename:0
140k_sample_test.4.mpileup.gz:67108864
140k_sample.5.list.samplename:0
140k_sample_test.5.mpileup.gz:67108864
...
```

编译 mpileup extractor，代码位于 BaseVar/MaxCompute/mpileup_extractor，编译后生成 mpileup-extractor-1.0.jar（依赖 java 和 maven）

```shell
mvn clean package
```

将所有 sample name 文件合并成一个，用于全局索引 sample name。

```bash
cat *.samplename | sort | uniq > all_names.txt
```

使用 odpscmd 完成外表的创建，注意 `'odps.properties.rolearn'='acs:ram::1817850323806830:role/aliyunodpsdefaultrole',` 这行根据实际情况会有不同（ [参考资料](https://help.aliyun.com/document_detail/45389.html)）。

```sql
-- 准备资源文件
add file all_names.txt;
add jar mpileup-extractor-1.0.jar;

set odps.task.major.version=2dot0_demo_flighting;
set odps.sql.planner.mode=lot;
set odps.sql.ddl.odps2=true;
set odps.sql.preparse.odps2=lot;

-- 创建外表， 注意 odps.properties.rolearn 的 value 根据实际情况会有不同
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

在 odpscmd 中验证外表创建成功：

```sql
odps@ huada_test>select * from oss_nifty limit 10;
+-------+-----+----------+-------------+----+----+----+
| chrid | pos | base_ref | sample_name | c1 | c2 | c3 |
+-------+-----+----------+-------------+----+----+----+
| chr11 | 5246595 | C        | jb28q3id11  | 1  | ,  | F  |
| chr11 | 5246595 | C        | fepqoeqjgr  | 0  | *  | *  |
| chr11 | 5246595 | C        | m864f7oqtg  | 0  | *  | *  |
| chr11 | 5246595 | C        | xbr7aszika  | 0  | *  | *  |
| chr11 | 5246595 | C        | x8qacbh61i  | 0  | *  | *  |
| chr11 | 5246595 | C        | qav498y893  | 0  | *  | *  |
| chr11 | 5246595 | C        | 6g94b8bflp  | 0  | *  | *  |
| chr11 | 5246595 | C        | f7pyq583gi  | 0  | *  | *  |
| chr11 | 5246595 | C        | ik424g3lsk  | 0  | *  | *  |
| chr11 | 5246595 | C        | qi8ibe1kwq  | 0  | *  | *  |
+-------+-----+----------+-------------+----+----+----+

```

### 构造百万列宽表（MaxCompute 内置表）

```sql
-- 创建 Python UDF 用于宽表聚合
add py expand_udf.py;
create function expand as expand_udf.Expand using expand_udf.py,all_names.txt;    

-- 创建宽表
create table nifty_expand as select chrid, pos, base_ref, expand(sample_name, c1, c2, c3) as one from nifty group by chrid, pos, base_ref;
```

SQL 的 group by 字句将所有相同点位的数据汇聚到一起，交由 expand_udf.py 定义的 expand 函数（[参考 MaxCompute 的聚合函数](https://help.aliyun.com/document_detail/27867.html)），函数中按全局索引 all_names.txt 指示的位置进行填充。

这样聚合出来的宽表具有以下特征：

- 只有原始数据覆盖到的点位才有记录
- 某个点位的某个 sample 如果没有数据，以 `0 * *` 填充

```sql
+------------+------------+------------+------------+
| chrid      | pos        | base_ref   | one        |
+------------+------------+------------+------------+
| chr11      | 5246595    | C          | 0      *       *       0       *       *       0       *       *	（后略）
+------------+------------+------------+------------+ 
```

## 基于 MaxCompute 的 BaseVar 算法

经过上述步骤，输入数据已经被整理为一个点位 100w sample 在一行的形式。因此基于 MaxCompute 的 BaseVar 算法实现相对简单，去掉了原单机代码中 adhoc 从各个地方随机读入数据并对齐的代码，精简了 basevar.py、basetype.py、mpileup.py 和 algorithm.py，表现为一个 MaxCompute UDF，并通过 UDF 的第一个参数（vcf 或 coverage）来指示 UDF 运行于 coverage 输出模式或 vcf 输出模式。

### 在 MaxCompute 中运行

因为 BaseVar 算法依赖了 scipy，而 scipy 默认并不在 MaxCompute 的运行环境当中，因此需要将 scipy 的包也当作 UDF 所需资源上传至 MaxCompute：

```shell
wget http://mirrors.aliyun.com/pypi/packages/ae/94/28ca6f9311e2351bb68da41ff8c1bc8f82bb82791f2ecd34efa953e60576/scipy-0.19.0-cp27-cp27m-manylinux1_x86_64.whl#md5=0e49f7fc8d31c1c79f0a4d63b29e8a1f
```

在 odpscmd 中创建资源及 UDF：

```sql
add archive scipy-0.19.1-cp27-cp27m-manylinux1_x86_64.whl as scipy.zip;
add py basevar.py;
add py basetype.py;
add py mpileup.py;
add py algorithm.py;
create function basevar as basevar.BaseVar using mpileup.py,basevar.py,basetype.py,algorithm.py,scipy.zip;
```

执行 SQL，产生 coverage 表

```sql
-- 创建 coverage 表
create table nifty_cvg (line string);

-- 执行 SQL 生成 coverage 表数据
set odps.task.major.version=2dot0_demo_flighting;
set odps.sql.planner.mode=lot;
set odps.sql.ddl.odps2=true;
set odps.sql.preparse.odps2=lot;

set odps.pypy.enabled=false;
set odps.isolation.session.enable = true;

insert overwrite table nifty_cvg select basevar('coverage', chrid, pos, base_ref, one) as cvg from nifty_expand;
```

执行 SQL，产生 vcf 表

```sql
-- 创建 vcf 表
create table nifty_vcf (line string);

-- 执行 SQL 生成 coverage 表数据
set odps.task.major.version=2dot0_demo_flighting;
set odps.sql.planner.mode=lot;
set odps.sql.ddl.odps2=true;
set odps.sql.preparse.odps2=lot;

set odps.pypy.enabled=false;
set odps.isolation.session.enable = true;

insert overwrite table nifty_vcf select basevar('vcf', chrid, pos, base_ref, one) as vcf from nifty_expand;
```

### 单机运行

为了 debug 和性能调优方便，特在 BaseVar UDF 中也添加了单机运行模式。单机运行模式下，可以直接将 MaxCompute 产生的宽表下载为本地文件作为程序的输入：

```sql
tunnel download nifty_expand nifty_expand.txt;
```

本地运行：

```shell
# vcf 模式
./basevar.py vcf nifty_expand.txt
# coverage 模式
./basevar.py vcf nifty_expand.txt
```

性能调优：

```shell
python -m cProfile -o perf.out ./basevar.py vcf nifty_expand.txt
```

对于产生的 perf.out，在 python 交互环境中：

```python
>>> import pstats
>>> p = pstats.Stats('perf.out')
>>> p.sort_stats('cumulative').print_stats(100)
Tue Jul 18 17:54:00 2017    perf.out

         76145918 function calls (76136563 primitive calls) in 187.993 seconds

   Ordered by: cumulative time
   List reduced from 6005 to 100 due to restriction <100>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    5.177    5.177  188.048  188.048 ./basevar.py:4(<module>)
      714   54.814    0.077  181.379    0.254 ./basevar.py:29(process)
      714   45.833    0.064   58.693    0.082 ./basetype.py:23(__init__)
  5776768   13.536    0.000   42.325    0.000 ./mpileup.py:173(first_base)
   643743   13.129    0.000   13.129    0.000 {numpy.core.multiarray.array}
       38    2.720    0.072   12.334    0.325 ./basevar.py:122(_out_vcf_line)
 11554035   11.050    0.000   11.051    0.000 {method 'sub' of '_sre.SRE_Pattern' objects}
  5776768    5.852    0.000   10.394    0.000 ./mpileup.py:154(scan_indel)
  5776768    3.127    0.000    9.407    0.000 ./mpileup.py:15(rmStartEnd)
       38    0.660    0.017    8.306    0.219 ./algorithm.py:85(strand_bias)
     3739    7.997    0.002    7.997    0.002 {method 'split' of 'str' objects}
       38    0.005    0.000    7.605    0.200 /home/tops/lib/python2.7/site-packages/scipy/stats/stats.py:3036(fisher_exact)
       32    0.280    0.009    7.534    0.235 /home/tops/lib/python2.7/site-packages/scipy/stats/stats.py:3117(binary_search)
  5776768    2.640    0.000    7.406    0.000 ./mpileup.py:32(rmIndel)
    85496    2.253    0.000    7.285    0.000 /home/tops/lib/python2.7/site-packages/scipy/stats/_distn_infrastructure.py:2812(pmf)
  5777114    4.544    0.000    4.544    0.000 {method 'search' of '_sre.SRE_Pattern' objects}
```

