package com.aliyun.odps.exttable.bam;

import com.aliyun.odps.udf.Extractor;
import com.aliyun.odps.udf.OdpsStorageHandler;
import com.aliyun.odps.udf.Outputer;

public class BamHandler extends OdpsStorageHandler {

    @Override
    public Class<? extends Extractor> getExtractorClass() {
        return BamExtractor.class;
    }

    @Override
    public Class<? extends Outputer> getOutputerClass() {
        return null;
    }
}
