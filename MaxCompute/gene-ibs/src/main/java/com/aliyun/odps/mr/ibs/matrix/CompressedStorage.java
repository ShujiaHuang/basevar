package com.aliyun.odps.mr.ibs.matrix;

public class CompressedStorage {

    public final int nSample, nPos;
    public final int[] samplePointers, posIndexes;
    public final double[] data;

    public CompressedStorage(int nSample, int nPos, MatrixElement[] elems) {
        this.nSample = nSample;
        this.nPos = nPos;
        samplePointers = new int[nSample + 1];
        posIndexes = new int[elems.length];
        data = new double[elems.length];

        MatrixElement prev = null;
        int curSample = 0;
        for (int i = 0; i < elems.length; i++) {
            MatrixElement e = elems[i];
            if (prev != null) {
                if (prev.sample > e.sample || (prev.sample == e.sample && prev.pos >= e.pos)) {
                    throw new IllegalArgumentException("Matrix elements are not in order (sorted by sample, then pos");
                }
            }
            while (e.sample != curSample) {
                curSample++;
                samplePointers[curSample] = i;
            }
            posIndexes[i] = e.pos;
            data[i] = e.value;
            prev = e;
        }
        curSample++;
        while (curSample <= nSample) {
            samplePointers[curSample] = elems.length;
            curSample++;
        }
    }
}
