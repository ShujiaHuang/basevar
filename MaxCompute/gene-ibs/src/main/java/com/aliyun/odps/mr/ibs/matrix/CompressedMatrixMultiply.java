package com.aliyun.odps.mr.ibs.matrix;

import java.util.ArrayList;
import java.util.List;

public class CompressedMatrixMultiply {
    public static MatrixElement[] multiply(CompressedStorage a, CompressedStorage b, double[] interceptionFm,
        double threshold) {
        List<int[]> outputRows = new ArrayList<>();
        List<int[]> outputCols = new ArrayList<>();
        List<double[]> outputVals = new ArrayList<>();
        List<int[]> outputCount = new ArrayList<>();
        int chunkSize = 1024 * 128;
        int[] curOutputRows = new int[chunkSize];
        int[] curOutputCols = new int[chunkSize];
        double[] curOutputVals = new double[chunkSize];
        int[] curOutputCount = new int[chunkSize];

        int nRows = a.nSample;
        int nCols = b.nSample;
        int[] aSamplePointers = a.samplePointers;
        int[] bSamplePointers = b.samplePointers;
        int[] aPosIndexes = a.posIndexes;
        int[] bPosIndexes = b.posIndexes;
        double[] aData = a.data;
        double[] bData = b.data;
        int i = 0;
        for (int r = 0; r < nRows; r++) {
            for (int c = 0; c < nCols; c++) {
                double sum = 0.0;
                int count = 0;
                int i1 = aSamplePointers[r];
                int i2 = bSamplePointers[c];
                int end1 = aSamplePointers[r+1];
                int end2 = bSamplePointers[c+1];
                while (i1 < end1 && i2 < end2) {
                    if (aPosIndexes[i1] == bPosIndexes[i2]) {
                        double fm = interceptionFm[aPosIndexes[i1]];
                        if (fm > threshold && fm < 1 - threshold) {
                            double tmp = (aData[i1] - fm) * (bData[i2] - fm) / fm / (1 - fm);
                            if (!Double.isNaN(tmp)) {
                                sum += tmp;
                                count++;
                            }
                        }
                        i1++; i2++;
                    } else if (aPosIndexes[i1] < bPosIndexes[i2]) {
                        i1++;
                    } else {
                        i2++;
                    }
                }
                if (count != 0) {
                    curOutputRows[i] = r;
                    curOutputCols[i] = c;
                    curOutputVals[i] = sum;
                    curOutputCount[i] = count;
                    i++;
                    if (i == chunkSize) {
                        outputRows.add(curOutputRows);
                        outputCols.add(curOutputCols);
                        outputVals.add(curOutputVals);
                        outputCount.add(curOutputCount);
                        curOutputRows = new int[chunkSize];
                        curOutputCols = new int[chunkSize];
                        curOutputVals = new double[chunkSize];
                        curOutputCount = new int[chunkSize];
                        i = 0;
                    }
                }
            }
        }
        MatrixElement[] result = new MatrixElement[outputRows.size() * chunkSize + i];
        int outI = 0;
        for (int j = 0; j < outputRows.size(); j++) {
            int[] rows = outputRows.get(j);
            int[] cols = outputCols.get(j);
            double[] vals = outputVals.get(j);
            int[] counts = outputCount.get(j);
            for (int k = 0; k < chunkSize; k++) {
                int row = rows[k];
                int col = cols[k];
                double val = vals[k];
                int count = counts[k];
                result[outI++] = new MatrixElement(row, col, val, count);
            }
        }
        for (int k = 0; k < i; k++) {
            int row = curOutputRows[k];
            int col = curOutputCols[k];
            double val = curOutputVals[k];
            int count = curOutputCount[k];
            result[outI++] = new MatrixElement(row, col, val, count);
        }

        return result;
    }
}
