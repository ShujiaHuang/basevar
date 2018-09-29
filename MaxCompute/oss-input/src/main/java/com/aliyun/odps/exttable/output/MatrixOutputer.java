package com.aliyun.odps.exttable.output;

import com.aliyun.odps.data.Record;
import com.aliyun.odps.io.OutputStreamSet;
import com.aliyun.odps.udf.DataAttributes;
import com.aliyun.odps.udf.ExecutionContext;
import com.aliyun.odps.udf.Outputer;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.Charset;
import java.util.zip.GZIPOutputStream;

public class MatrixOutputer extends Outputer {
    private long currentI = -1l;
    private Writer writer = null;
    private OutputStreamSet outputStreamSet = null;

    public void setup(ExecutionContext executionContext, OutputStreamSet outputStreamSet, DataAttributes dataAttributes) {
        this.outputStreamSet = outputStreamSet;
    }

    public void output(Record record) throws IOException {
        if (writer == null) {
            long partitionKey = record.getBigint(0);
            writer = new OutputStreamWriter(new GZIPOutputStream(outputStreamSet.next(partitionKey + ".gz")), "UTF-8");
        }
        long i = record.getBigint(1);
        double val = record.getDouble(3);
        if (currentI != i) {
            if (currentI != -1) {
                writer.write("\n");
            }
            currentI = i;
        }
        writer.write(String.format("%.6f\t", val));
    }

    public void close() throws IOException {
        if (writer != null) {
            writer.close();
        }
    }
}
