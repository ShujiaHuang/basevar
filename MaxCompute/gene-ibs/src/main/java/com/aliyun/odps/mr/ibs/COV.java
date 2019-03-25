package com.aliyun.odps.mr.ibs;

import com.aliyun.odps.OdpsException;
import com.aliyun.odps.data.TableInfo;
import com.aliyun.odps.mapred.JobClient;
import com.aliyun.odps.mapred.conf.JobConf;
import com.aliyun.odps.mapred.utils.InputUtils;
import com.aliyun.odps.mapred.utils.OutputUtils;
import com.aliyun.odps.mapred.utils.SchemaUtils;
import com.aliyun.odps.mr.ibs.cov.mapred.CovMapper;
import com.aliyun.odps.mr.ibs.cov.mapred.CovPartitioner;
import com.aliyun.odps.mr.ibs.cov.mapred.CovReducer;

import java.io.IOException;

import static java.lang.System.exit;

/**
 * Created by tianli on 2017/8/8.
 */
public class COV {
    public static void main(String[] args) throws IOException, OdpsException {
        if (args.length != 7) {
            System.err.println("args: input_table output_table sample_block_count sample_count pos_block_count pos_count fm_threshold");
            exit(1);
        }

        String inputTable = args[0];
        String outputTable = args[1];
        int sampleBlockCount = Integer.valueOf(args[2]);
        int sampleCount = Integer.valueOf(args[3]);
        int posBlockCount = Integer.valueOf(args[4]);
        int posCount = Integer.valueOf(args[5]);
        String fmThreshold = args[6];

        JobConf jobConf = new JobConf();
        jobConf.setInt("sample_block_count", sampleBlockCount);
        jobConf.setInt("sample_count", sampleCount);
        jobConf.setInt("pos_block_count", posBlockCount);
        jobConf.setInt("pos_count", posCount);
        jobConf.set("fm_threshold", fmThreshold);

        jobConf.setMapperClass(CovMapper.class);
        jobConf.setReducerClass(CovReducer.class);
        jobConf.setPartitionerClass(CovPartitioner.class);
        jobConf.setMapOutputKeySchema(SchemaUtils.fromString(
            "iblock:bigint,jblock:bigint,kblock:bigint,m:bigint,sample:bigint,pos:bigint"));
        jobConf.setMapOutputValueSchema(SchemaUtils.fromString("val:double"));
        jobConf.setOutputKeySortColumns(new String[]{ "iblock", "jblock", "kblock" });
        jobConf.setOutputGroupingColumns(new String[]{ "iblock", "jblock", "kblock" });

        jobConf.setNumReduceTasks((sampleBlockCount + 1) * sampleBlockCount / 2);

        InputUtils.addTable(TableInfo.builder().tableName(inputTable).build(), jobConf);
        OutputUtils.addTable(TableInfo.builder().tableName(outputTable).build(), jobConf);

        JobClient.runJob(jobConf);
    }

}
