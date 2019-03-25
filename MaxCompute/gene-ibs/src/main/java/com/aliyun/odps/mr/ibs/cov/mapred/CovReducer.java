package com.aliyun.odps.mr.ibs.cov.mapred;

import com.aliyun.odps.data.Record;
import com.aliyun.odps.mapred.ReducerBase;
import com.aliyun.odps.mapred.conf.JobConf;
import com.aliyun.odps.mr.ibs.matrix.CompressedMatrixMultiply;
import com.aliyun.odps.mr.ibs.matrix.CompressedStorage;
import com.aliyun.odps.mr.ibs.matrix.MatrixBuilder;
import com.aliyun.odps.mr.ibs.matrix.MatrixElement;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class CovReducer extends ReducerBase {

    public final static String SAMPLE = "sample";
    public final static String POS = "pos";
    public final static String VAL = "val";

    private Record result;
    private int sib = -1, sjb = -1, skb = -1;
    private int sampleBlockSize;

    private List<Integer> sampleIdxA;
    private List<List<Integer>> posIdxA;
    private List<List<Double>> dataA;
    private double[] acc;
    private int[] count;
    private double[] fm;
    private double fmThreshold;

    private long timeFm = 0l, timeM0 = 0l, timeM1 = 0l, timeWriteBlock = 0l, timeCalculate = 0l;
    private long accFm = 0l, accM0 = 0l, accM1 = 0l, accWriteBlock = 0l, accCalculate = 0l;
    private long readAcc = 0l, filterAcc = 0l, compressAcc = 0l, multiplyAcc = 0l, updateAcc = 0l;

    private static String getMillTimeStr(long value) {
        return TimeUnit.NANOSECONDS.toMillis(value) + "ms";
    }

    private static String getMicroTimeStr(long value) {
        return TimeUnit.NANOSECONDS.toMicros(value) + "us";
    }


    @Override
    public void setup(TaskContext context) throws IOException {
        result = context.createOutputRecord();

        JobConf jobConf = context.getJobConf();
        int sampleBlockCount = jobConf.getInt("sample_block_count", 60);
        int sampleCount = jobConf.getInt("sample_count", 1000000);
        int posCount = jobConf.getInt("pos_count", 6000000);
        fmThreshold = Double.valueOf(jobConf.get("fm_threshold", "0.01"));
        System.out.println("fmThreshold: " + fmThreshold);
        sampleBlockSize = (sampleCount - 1) / sampleBlockCount + 1;
        int resSize = sampleBlockSize * sampleBlockSize;
        acc = new double[resSize];
        count = new int[resSize];
        fm = new double[posCount];
    }

    private int getResIdx(int i, int j) {
        return i * sampleBlockSize + j;
    }

    public void processBlock(Record key, Iterator<Record> values) {
        long start = System.nanoTime();
        List<MatrixElement> leftElems = new ArrayList<>();
        List<MatrixElement> rightElems = new ArrayList<>();
        readElements(leftElems, rightElems, key, values); // This method also updates fm array
        long readComplete = System.nanoTime();
        Set<Integer> intersection = leftElems.stream().map(e -> e.pos).collect(Collectors.toSet());
        intersection.retainAll(rightElems.stream().map(e -> e.pos).collect(Collectors.toSet()));
        List<MatrixElement> leftIncludedElements = leftElems.stream().filter(e -> intersection.contains(e.pos)).collect(Collectors.toList());
        List<MatrixElement> rightIncludedElements = rightElems.stream().filter(e -> intersection.contains(e.pos)).collect(Collectors.toList());

        List<Integer> leftSamples = leftIncludedElements.stream().map(e -> e.sample).distinct().collect(Collectors.toList());
        List<Integer> intersectionPos = intersection.stream().collect(Collectors.toList());
        List<Integer> rightSamples = rightIncludedElements.stream().map(e -> e.sample).distinct().collect(Collectors.toList());

        Map<Integer, Integer> leftSampleMap =
            IntStream.range(0, leftSamples.size())
                .boxed()
                .collect(Collectors.toMap(i -> leftSamples.get(i), i -> i));
        Map<Integer, Integer> posMap =
            IntStream.range(0, intersectionPos.size())
                .boxed()
                .collect(Collectors.toMap(i -> intersectionPos.get(i), i -> i));
        Map<Integer, Integer> rightSampleMap =
            IntStream.range(0, rightSamples.size())
                .boxed()
                .collect(Collectors.toMap(i -> rightSamples.get(i), i -> i));
        double[] interceptionFm = intersectionPos.stream().mapToDouble(i -> fm[i]).toArray();

        long sampleFilterComplete = System.nanoTime();
        CompressedStorage leftMatrix =
            new MatrixBuilder(leftIncludedElements.stream().map(
                e -> new MatrixElement(leftSampleMap.get(e.sample), posMap.get(e.pos), e.value)
            ).collect(Collectors.toList())).buildCompressedStorage();
        CompressedStorage rightMatrix =
            new MatrixBuilder(rightIncludedElements.stream().map(
                e -> new MatrixElement(rightSampleMap.get(e.sample), posMap.get(e.pos), e.value)
            ).collect(Collectors.toList())).buildCompressedStorage();
        statisticMatrix(leftMatrix.samplePointers);
        statisticMatrix(rightMatrix.samplePointers);
        long buildMatrixComplete = System.nanoTime();
        MatrixElement[] result = CompressedMatrixMultiply.multiply(leftMatrix, rightMatrix, interceptionFm, fmThreshold);
        long multiplyComplete = System.nanoTime();
        updateCache(leftSamples, rightSamples, result);
        long updateCacheComplete = System.nanoTime();

        logTime(start, readComplete, sampleFilterComplete, buildMatrixComplete, multiplyComplete, updateCacheComplete);
    }

    private void statisticMatrix(int[] samplePointers) {
        int max = 0;
        int min = Integer.MAX_VALUE;
        double avg = samplePointers[samplePointers.length - 1] / samplePointers.length;
        for (int i = 0; i < samplePointers.length - 1; i++) {
            int tmp = samplePointers[i+1] - samplePointers[i];
            if (tmp > max) { max = tmp; }
            if (tmp < min) { min = tmp; }
        }
        System.out.println(String.format("max: %s, min: %s, avg: %s", max, min, avg));
    }

    private void logTime(long start, long read, long filter, long compress, long multiply, long update) {
        long r = read - start;
        long f = filter - read;
        long c = compress - filter;
        long m = multiply - compress;
        long u = update - multiply;
        System.out.println(String.format("Read-%s, filter-%s, compress-%s, mult-%s, update-%s",
            getMillTimeStr(r),
            getMillTimeStr(f),
            getMillTimeStr(c),
            getMillTimeStr(m),
            getMillTimeStr(u))
        );
        readAcc += r;
        filterAcc += f;
        compressAcc += c;
        multiplyAcc += m;
        updateAcc += u;
    }

    private void logAccTime() {
        System.out.println(String.format("ACC: Read-%s, filter-%s, compress-%s, mult-%s, update-%s",
            getMillTimeStr(readAcc),
            getMillTimeStr(filterAcc),
            getMillTimeStr(compressAcc),
            getMillTimeStr(multiplyAcc),
            getMillTimeStr(updateAcc))
        );
    }

    private void updateCache(List<Integer> leftSamples, List<Integer> rightSamples, MatrixElement[] result) {
        for (MatrixElement element : result) {
            int curCount = element.count;
            if (curCount <= 0) {
                continue;
            }
            int row = leftSamples.get(element.sample); // sample stores row
            int col = rightSamples.get(element.pos); // pos stores col
            double sum = element.value;
            int resIdx = getResIdx(row, col);
            acc[resIdx] += sum;
            count[resIdx] += curCount;
        }
    }

    private void readElements(List<MatrixElement> leftElems, List<MatrixElement> rightElems, Record key,
        Iterator<Record> values) {

        while (values.hasNext()) {
            Record value = values.next();
            int m = key.getBigint("m").intValue();
            MatrixElement elem = new MatrixElement(key.getBigint(SAMPLE).intValue(), key.getBigint(POS).intValue(),
                value.getDouble(VAL));
            if (m == 0) {
                leftElems.add(elem);
            } else if (m == 1) {
                rightElems.add(elem);
            } else if (m == -1) {
                int pos = key.getBigint(POS).intValue();
                fm[pos] = value.getDouble(VAL);
            } else {
                throw new RuntimeException("Unexpected m: " + m);
            }
        }
    }

    @Override
    public void reduce(Record key, Iterator<Record> values, TaskContext context) throws  IOException {
        int ib = key.getBigint("iblock").intValue();
        int jb = key.getBigint("jblock").intValue();
        int kb = key.getBigint("kblock").intValue();
        if (ib != sib || jb != sjb) {
            if (sib != -1) {
                System.out.println(String.format("ib:%s,jb:%s,sib:%s,sjb:%s", ib, jb, sib, sjb));
                throw new RuntimeException("Unexpected change block.");
            }
            sib = ib;
            sjb = jb;
            skb = -1;
            // Zero C
        }
        System.out.println(String.format("[%s] kb: %s" , new Date().toLocaleString(), kb));
        processBlock(key, values);
    }

    @Override
    public void cleanup(TaskContext context) throws IOException {
        writeBlockResult(context);
        logAccTime();
        System.out.println(String.format("Fm: %s, M0: %s, M1: %s, Cal: %s, Write: %s",
            getMillTimeStr(accFm),
            getMillTimeStr(accM0),
            getMillTimeStr(accM1),
            getMillTimeStr(accCalculate),
            getMillTimeStr(accWriteBlock)));
    }

    private void writeBlockResult(TaskContext context) throws IOException {
        long start = System.nanoTime();
        int ibase = sib * sampleBlockSize;
        int jbase = sjb * sampleBlockSize;
        for (int i = 0; i < sampleBlockSize; i++) {
            for (int j = 0; j < sampleBlockSize; j++) {
                // Emit C[i,j]
                int resIdx = getResIdx(i, j);
                if (count[resIdx] == 0) {
                    continue;
                }
                int realI = i + ibase;
                int realJ = j + jbase;
                if (realI > realJ) { continue; }
                result.setBigint(0, (long) realI);
                result.setBigint(1, (long) realJ);
                result.setDouble(2, acc[resIdx] / count[resIdx]);
                result.setBigint(3, (long) count[resIdx]);
                context.write(result);
                count[resIdx] = 0;
                acc[resIdx] = 0;
            }
        }
        timeWriteBlock = System.nanoTime() - start;
    }
}
