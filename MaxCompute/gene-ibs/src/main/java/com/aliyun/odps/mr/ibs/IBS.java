package com.aliyun.odps.mr.ibs;

import com.aliyun.odps.OdpsException;
import com.aliyun.odps.data.Record;
import com.aliyun.odps.data.TableInfo;
import com.aliyun.odps.mapred.*;
import com.aliyun.odps.mapred.conf.JobConf;
import com.aliyun.odps.mapred.utils.InputUtils;
import com.aliyun.odps.mapred.utils.OutputUtils;
import com.aliyun.odps.mapred.utils.SchemaUtils;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;
import java.util.List;

import static java.lang.System.exit;

/**
 * Created by tianli on 2017/8/8.
 */
public class IBS {

    public static class IBSMapper extends MapperBase {
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
            String val = record.getString("val");

            value.setString(0, val);

            long iblock = sample / sampleBlockSize;
            key.setBigint("iblock", iblock);
            key.setBigint("kblock", pos / posBlockSize);
            key.setBigint("m", 0l);
            key.setBigint("sample", sample % sampleBlockSize);
            key.setBigint("pos", pos % posBlockSize);
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

    public static class IBSPartitioner extends Partitioner {

        private int sampleBlockCount;

        @Override
        public void configure(JobConf jobConf) {
            sampleBlockCount = jobConf.getInt("sample_block_count", 60);
        }

        @Override
        public int getPartition(Record key, Record value, int numPartitions) {
            int ib = key.getBigint("iblock").intValue();
            int jb = key.getBigint("jblock").intValue();
            return (ib * sampleBlockCount + jb) % numPartitions; // Modify, maybe doc is wrong, need test
        }
    }

    public static class IBSReducer extends ReducerBase {
        private Record result;
        private int sib = -1, sjb = -1, skb = -1;
        private int sampleBlockSize;

        private List<Integer> sampleIdxA;
        private List<List<Integer>> posIdxA;
        private List<List<String>> dataA;
        private int[] total;
        private int[] same;

        @Override
        public void setup(TaskContext context) throws IOException {
            result = context.createOutputRecord();

            JobConf jobConf = context.getJobConf();
            int sampleBlockCount = jobConf.getInt("sample_block_count", 60);
            int sampleCount = jobConf.getInt("sample_count", 1000000);
            sampleBlockSize = (sampleCount - 1) / sampleBlockCount + 1;
            int resSize = sampleBlockSize * sampleBlockSize;
            total = new int[resSize];
            same = new int[resSize];
        }

        private int getResIdx(int i, int j) {
            return i * sampleBlockSize + j;
        }



        @Override
        public void reduce(Record key, Iterator<Record> values, TaskContext context) throws  IOException {
            int ib = key.getBigint("iblock").intValue();
            int jb = key.getBigint("jblock").intValue();
            int kb = key.getBigint("kblock").intValue();
            int m = key.getBigint("m").intValue();
            System.err.println(
                String.format("[%s] ib=%s, jb=%s, kb=%s, m=%s", new Date().toLocaleString(), ib, jb, kb, m));
            if (ib != sib || jb != sjb) {
                if (sib != -1) {
                    writeBlockResult(context);
                }
                sib = ib;
                sjb = jb;
                skb = -1;
                // Zero C
            }
            if (m == 0) {
                skb = kb;
                dataA = new ArrayList<List<String>>();
                sampleIdxA = new ArrayList<Integer>();
                posIdxA = new ArrayList<List<Integer>>();
                List<String> currentDataList = null;
                List<Integer> currentPosList = null;
                int currentSample = -1;

                while (values.hasNext()) {
                    Record value = values.next();
                    int sample = key.getBigint("sample").intValue();
                    int pos = key.getBigint("pos").intValue();
                    String c = value.getString("val");
                    if (currentSample != sample) {
                        currentSample = sample;

                        currentDataList = new ArrayList<String>();
                        dataA.add(currentDataList);
                        currentPosList = new ArrayList<Integer>();
                        posIdxA.add(currentPosList);

                        sampleIdxA.add(sample);
                    }
                    currentDataList.add(c);
                    currentPosList.add(pos);
                }
            }
            if (m == 1) {
                if (kb != skb) { return; }

                List<String> currentBDataList = null;
                List<Integer> currentBPosList = null;
                int currentBSample = -1;
                while (values.hasNext()) {
                    Record value = values.next();
                    int sample = key.getBigint("sample").intValue();
                    int pos = key.getBigint("pos").intValue();
                    String c = value.getString("val");
                    if (currentBSample != sample) {
                        if (currentBSample != -1) {
                            calculateSample(total, same, currentBSample, currentBDataList, currentBPosList);
                        }
                        currentBDataList = new ArrayList<String>();
                        currentBPosList = new ArrayList<Integer>();
                        currentBSample = sample;
                    }
                    currentBDataList.add(c);
                    currentBPosList.add(pos);
                }
                calculateSample(total, same, currentBSample, currentBDataList, currentBPosList);
            }
        }

        @Override
        public void cleanup(TaskContext context) throws IOException {
            writeBlockResult(context);
        }

        private void writeBlockResult(TaskContext context) throws IOException {
            int ibase = sib * sampleBlockSize;
            int jbase = sjb * sampleBlockSize;
            for (int i = 0; i < sampleBlockSize; i++) {
                for (int j = 0; j < sampleBlockSize; j++) {
                    // Emit C[i,j]
                    int resIdx = getResIdx(i, j);
                    if (total[resIdx] == 0) {
                        continue;
                    }
                    int realI = i + ibase;
                    int realJ = j + jbase;
                    if (realI > realJ) { continue; }
                    result.setBigint(0, (long) realI);
                    result.setBigint(1, (long) realJ);
                    result.setDouble(2, 1 - 1.0 * same[resIdx] / total[resIdx]);
                    result.setBigint(3, (long) total[resIdx]);
                    context.write(result);
                    same[resIdx] = 0;
                    total[resIdx] = 0;
                }
            }
        }

        private void calculateSample(int[] total, int[] same, int sampleB, List<String> listB, List<Integer> posB) {
            for (int i = 0; i < sampleIdxA.size(); i++) {
                int sampleA = sampleIdxA.get(i);
                int resIdx = getResIdx(sampleA, sampleB);
                List<String> listA = dataA.get(i);
                List<Integer> posA = posIdxA.get(i);

                int a = 0, b = 0;
                int sizeA = listA.size();
                int sizeB = listB.size();
                while (a < sizeA && b < sizeB) {
                    int currentPosA = posA.get(a);
                    int currentPosB = posB.get(b);
                    if (currentPosA == currentPosB) {
                        ++total[resIdx];
                        String aVal = listA.get(a);
                        String bval = listB.get(b);
                        if (aVal.equalsIgnoreCase(bval)) {
                            ++same[resIdx];
                        }
                        ++a; ++b;
                        continue;
                    }
                    if (currentPosA < currentPosB) {
                        ++a;
                        continue;
                    }
                    if (currentPosA > currentPosB) {
                        ++b;
                    }
                }
            }
        }
    }

    public static void main(String[] args) throws IOException, OdpsException {
        if (args.length != 6) {
            System.err.println("args: input_table output_table sample_block_count sample_count pos_block_count pos_count");
            exit(1);
        }

        String inputTable = args[0];
        String outputTable = args[1];
        int sampleBlockCount = Integer.valueOf(args[2]);
        int sampleCount = Integer.valueOf(args[3]);
        int posBlockCount = Integer.valueOf(args[4]);
        int posCount = Integer.valueOf(args[5]);

        JobConf jobConf = new JobConf();
        jobConf.setInt("sample_block_count", sampleBlockCount);
        jobConf.setInt("sample_count", sampleCount);
        jobConf.setInt("pos_block_count", posBlockCount);
        jobConf.setInt("pos_count", posCount);

        jobConf.setMapperClass(IBSMapper.class);
        jobConf.setReducerClass(IBSReducer.class);
        jobConf.setPartitionerClass(IBSPartitioner.class);
        jobConf.setMapOutputKeySchema(SchemaUtils.fromString(
            "iblock:bigint,jblock:bigint,kblock:bigint,m:bigint,sample:bigint,pos:bigint"));
        jobConf.setMapOutputValueSchema(SchemaUtils.fromString("val:string"));
        jobConf.setOutputKeySortColumns(new String[]{ "iblock", "jblock", "kblock", "m", "sample", "pos" });
        jobConf.setOutputGroupingColumns(new String[]{ "iblock", "jblock", "kblock", "m" });

        jobConf.setNumReduceTasks(sampleBlockCount * sampleBlockCount);

        InputUtils.addTable(TableInfo.builder().tableName(inputTable).build(), jobConf);
        OutputUtils.addTable(TableInfo.builder().tableName(outputTable).build(), jobConf);

        JobClient.runJob(jobConf);
    }

}
