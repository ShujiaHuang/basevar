package com.aliyun.odps.mr.ibs.matrix;

import com.aliyun.odps.Column;
import com.aliyun.odps.OdpsType;
import com.aliyun.odps.data.ArrayRecord;
import com.aliyun.odps.data.Record;
import com.aliyun.odps.mr.ibs.cov.mapred.CovReducer;

import java.util.Iterator;

public class TestCovReducer {
    public void testProcessBlock() {
        CovReducer r = new CovReducer();

        Column[] keyColumns = new Column[]{
            new Column("m", OdpsType.BIGINT),
            new Column("sample", OdpsType.BIGINT),
            new Column("pos", OdpsType.BIGINT)
        };
        Record key = new ArrayRecord(keyColumns);

        Column[] valueColumns = new Column[] {new Column("val", OdpsType.DOUBLE)};
        Record value = new ArrayRecord(valueColumns);

        double[] fm = new double[] {0.12, 0.13, 0.14, 0.15};
        String[] leftSamples = new String[] {"1,0,1,1", "0,1,0,0"};
        String[] rightSamples = new String[] { "1,1,0,1","0,1,0,1"};

        Iterator<Record> values = new Iterator<Record>() {
            int total = leftSamples.length + rightSamples.length + fm.length;
            int cur = 0, leftIdx = 0, rightIdx = 0, fmIdx = 0, posIdx = 0;
            String[] currentTokens = null;

            @Override
            public boolean hasNext() {
                return cur < 3;
            }

            @Override
            public Record next() {
                if (leftIdx < leftSamples.length) {
                    key.setBigint("m", 0l);
                    key.setBigint("sample", Long.valueOf(leftIdx));
                    key.set("pos", Long.valueOf(posIdx));

                    if (currentTokens == null) {
                        currentTokens = leftSamples[leftIdx].split(",");
                    }
                    value.set("val", currentTokens[posIdx]);
                    posIdx++;
                    if (posIdx == currentTokens.length) {
                        posIdx = 0;
                        leftIdx++;
                        if (leftIdx == leftSamples.length) {
                            cur++;
                        }
                    }
                } else if (rightIdx < rightSamples.length) {
                    key.setBigint("m", 1l);
                    key.setBigint("sample", Long.valueOf(rightIdx));
                }

                return value;
            }
        };

        r.processBlock(key, values);
    }
}
