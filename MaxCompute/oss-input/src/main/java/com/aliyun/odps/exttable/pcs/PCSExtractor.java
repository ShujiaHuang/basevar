package com.aliyun.odps.exttable.pcs;

import com.aliyun.odps.Column;
import com.aliyun.odps.OdpsType;
import com.aliyun.odps.data.ArrayRecord;
import com.aliyun.odps.data.Record;
import com.aliyun.odps.exttable.ibs.AvailableInputStream;
import com.aliyun.odps.io.InputStreamSet;
import com.aliyun.odps.io.SourceInputStream;
import com.aliyun.odps.udf.DataAttributes;
import com.aliyun.odps.udf.ExecutionContext;
import com.aliyun.odps.udf.Extractor;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class PCSExtractor extends Extractor {
    private InputStreamSet inputs;
    private ArrayRecord record;
    private String fileSuffix = null;
    private int sampleIdx = 0;
    private String sampleName = null;
    private boolean isFirstRead = true;
    private BufferedReader br = null;
    private String mode;

    @Override
    public void setup(ExecutionContext ctx, InputStreamSet inputs, DataAttributes attributes) {
        this.inputs = inputs;
        record = new ArrayRecord(getColumns());
        mode = attributes.getValueByKey("mode");
        if (mode == null) {
            mode = "PC";
        }
    }

    private Column[] getColumns() {
        Column[] columns = new Column[4];
        columns[0] = new Column("suffix", OdpsType.STRING);
        columns[1] = new Column("sample_idx", OdpsType.BIGINT);
        columns[2] = new Column("val", OdpsType.DOUBLE);
        columns[3] = new Column("sample_name", OdpsType.STRING);
        return columns;
    }

    private BufferedReader getNextReader() throws IOException {
        String fileName;
        SourceInputStream sourceInputStream;

        do {
            sourceInputStream = inputs.next();
            if (sourceInputStream == null) {
                return null;
            }

            fileName = sourceInputStream.getFileName();
            System.err.println(fileName);
        } while (!fileName.contains(mode));
        fileSuffix = fileName.substring(fileName.lastIndexOf('.') + 1);
        InputStream in = new AvailableInputStream(sourceInputStream);
        BufferedReader br = new BufferedReader(new InputStreamReader(in));
        br.readLine();
        sampleIdx = -1;
        return br;
    }

    private String[] getNextLine() throws IOException {
        String[] items;
        if (isFirstRead) {
            isFirstRead = false;
            br = getNextReader();
        }
        while (br != null) {
            String line = br.readLine();
            if (line == null) {
                br.close();
                br = getNextReader();
                continue;
            }
            items = line.split("\t");
            ++sampleIdx;
            return items;
        }
        return null;
    }

    private String getNextVal() throws IOException {
        String[] cols = getNextLine();
        if (cols == null) {
            return null;
        }
        if (mode.equalsIgnoreCase("PC")) {
            sampleName = cols[0];
            return cols[1];
        } else if (mode.equalsIgnoreCase("EIGEN")) {
            return cols[0];
        }
        throw new IllegalArgumentException("Unsupported mode: " + mode);
    }

    @Override
    public Record extract() throws IOException {
        String val = getNextVal();
        if (val == null) {
            return null;
        }
        record.setString("suffix", fileSuffix);
        record.setBigint("sample_idx", Long.valueOf(sampleIdx));
        record.setDouble("val", Double.valueOf(val));
        record.setString("sample_name", sampleName);
        return record;
    }

    @Override
    public void close() {

    }
}
