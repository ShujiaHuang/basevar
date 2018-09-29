package com.aliyun.odps.mr.ibs.matrix;

public class MatrixElement {
    public final int sample;
    public final int pos;
    public final double value;
    public final int count;

    public MatrixElement(int sample, int pos, double value) {
        this(sample, pos, value, -1);
    }

    public MatrixElement(int sample, int pos, double value, int count) {
        this.sample = sample;
        this.pos = pos;
        this.value = value;
        this.count = count;
    }

    public int getSample() {return this.sample;}
    public int getPos() {return this.pos;}

    @Override
    public String toString() {
        return "MatrixElement(" + sample + "," + pos + "," + value + ")";
    }

}
