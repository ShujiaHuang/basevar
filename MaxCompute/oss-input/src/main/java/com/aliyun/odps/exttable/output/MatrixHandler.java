package com.aliyun.odps.exttable.output;

import com.aliyun.odps.udf.Extractor;
import com.aliyun.odps.udf.OdpsStorageHandler;
import com.aliyun.odps.udf.Outputer;

public class MatrixHandler extends OdpsStorageHandler {
    public Class<? extends Extractor> getExtractorClass() {
        return null;
    }

    public Class<? extends Outputer> getOutputerClass() {
        return MatrixOutputer.class;
    }
}
