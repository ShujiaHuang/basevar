package com.aliyun.odps.exttable.pcs;

import com.aliyun.odps.udf.Extractor;
import com.aliyun.odps.udf.OdpsStorageHandler;
import com.aliyun.odps.udf.Outputer;

public class PCSHandler extends OdpsStorageHandler {
    @Override
    public Class<? extends Extractor> getExtractorClass() {
        return PCSExtractor.class;
    }

    @Override
    public Class<? extends Outputer> getOutputerClass() {
        return null;
    }
}
