package com.aliyun.odps.mr.ibs.matrix;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Stream;

public class MatrixBuilder {
    private List<MatrixElement> elements;
    public MatrixBuilder(List<MatrixElement> elements) {
        this.elements = elements;
    }

    private int getNSample() {
        return elements.size() == 0 ? 0 :
            elements.stream().map(e -> e.sample).max(Comparator.comparingInt(i -> i)).get() + 1;
    }

    private int getNPos() {
        return elements.size() == 0 ? 0 :
            elements.stream().map(e -> e.pos).max(Comparator.comparingInt(i -> i)).get() + 1;
    }

    public CompressedStorage buildCompressedStorage() {
        MatrixElement[] arr =
            elements.stream().sorted(Comparator.comparing(MatrixElement::getSample).thenComparing(MatrixElement::getPos))
                .toArray(MatrixElement[]::new);
        return new CompressedStorage(getNSample(), getNPos(), arr);
    }
}
