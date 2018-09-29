package com.aliyun.odps.mr.ibs.cov.mapred;

import com.aliyun.odps.data.Record;
import com.aliyun.odps.mapred.Partitioner;
import com.aliyun.odps.mapred.conf.JobConf;

public class CovPartitioner extends Partitioner {
    private int sampleBlockCount;

    @Override
    public void configure(JobConf jobConf) {
        sampleBlockCount = jobConf.getInt("sample_block_count", 60);
    }

    @Override
    public int getPartition(Record key, Record value, int numPartitions) {
        int ib = key.getBigint("iblock").intValue();
        int jb = key.getBigint("jblock").intValue();
        return (ib * (2 * sampleBlockCount - ib - 1) / 2 + jb) % numPartitions;
        //return (ib * sampleBlockCount + jb) % numPartitions; // Modify, maybe doc is wrong, need test
    }
}
