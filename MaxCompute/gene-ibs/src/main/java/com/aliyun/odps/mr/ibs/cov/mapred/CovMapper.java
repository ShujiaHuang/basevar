package com.aliyun.odps.mr.ibs.cov.mapred;

import com.aliyun.odps.data.Record;
import com.aliyun.odps.mapred.MapperBase;
import com.aliyun.odps.mapred.conf.JobConf;

import java.io.IOException;

public class CovMapper extends MapperBase {
    private Record key, value;
    private int sampleBlockSize, posBlockSize, sampleBlockCount;

    @Override
    public void setup(TaskContext context) throws IOException {
        key = context.createMapOutputKeyRecord();
        value = context.createMapOutputValueRecord();

        JobConf jobConf = context.getJobConf();
        sampleBlockCount = jobConf.getInt("sample_block_count", 60); // Total block = blockCount ^ 2
        int sampleCount = jobConf.getInt("sample_count", 1000000);
        int posCount = jobConf.getInt("pos_count", 6000000);
        int posBlockCount = jobConf.getInt("pos_block_count", 60);

        sampleBlockSize = (sampleCount - 1) / sampleBlockCount + 1;
        posBlockSize = (posCount - 1) / posBlockCount + 1;
        System.out.println("TaskID: " + context.getTaskID().toString());
    }

    @Override
    public void map(long recordNum, Record record, TaskContext context) throws IOException {
        long sample = record.getBigint("sample_idx");
        long pos = record.getBigint("pos_idx");
        double val = record.getDouble("val");

        value.setDouble(0, val);

        key.setBigint("kblock", pos / posBlockSize);
        key.setBigint("sample", sample % sampleBlockSize);
        key.setBigint("pos", pos % posBlockSize);
        if (sample == -1) {
            key.setBigint("m", -1l);
            for (long i = 0; i < sampleBlockCount; i++) {
                for (long j = i; j < sampleBlockCount; j++) {
                    key.setBigint("iblock", i);
                    key.setBigint("jblock", j);
                    context.write(key, value);
                }
            }
            return;
        }

        long iblock = sample / sampleBlockSize;
        key.setBigint("iblock", iblock);
        key.setBigint("m", 0l);
        for (long j = 0; j < sampleBlockCount; j++) {
            if (iblock > j) { continue; }
            key.setBigint("jblock", j);
            context.write(key, value);
        }

        long jblock = sample / sampleBlockSize;
        key.set("jblock", jblock);
        key.setBigint("m", 1l);
        for (long i = 0; i < sampleBlockCount; i++) {
            if (i > jblock) { continue; }
            key.setBigint("iblock", i);
            context.write(key, value);
        }
    }
}
