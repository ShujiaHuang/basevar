package com.aliyun.odps.exttable.output;

import com.aliyun.odps.data.Record;
import com.aliyun.odps.io.OutputStreamSet;
import com.aliyun.odps.udf.DataAttributes;
import com.aliyun.odps.udf.ExecutionContext;
import com.aliyun.odps.udf.Outputer;

import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.zip.GZIPOutputStream;

public class BasevarOutputer extends Outputer {

    private OutputStreamSet outputStreamSet = null;
    private Writer writer = null;
    private int currentIdx = 0;

    @Override
    public void setup(ExecutionContext ctx, OutputStreamSet outputStreamSet, DataAttributes attributes) {
        this.outputStreamSet = outputStreamSet;
    }

    @Override
    public void output(Record record) throws IOException {
        if (writer == null) {
            long partitionKey = record.getBigint(0);
            writer = new OutputStreamWriter(new GZIPOutputStream(outputStreamSet.next(partitionKey + ".gz")), "UTF-8");
        }
        String line = record.getString(2);
        writer.write(line);
        writer.write("\n");
    }

    @Override
    public void close() throws IOException {
        if (writer != null) {
            writer.close();
        }
    }
}
