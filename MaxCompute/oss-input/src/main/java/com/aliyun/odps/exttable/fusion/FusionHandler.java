package com.aliyun.odps.exttable.fusion;

import com.aliyun.odps.udf.Extractor;
import com.aliyun.odps.udf.OdpsStorageHandler;
import com.aliyun.odps.udf.Outputer;

public class FusionHandler extends OdpsStorageHandler {

    @Override
    public Class<? extends Extractor> getExtractorClass() {
        return FusionExtractor.class;
    }

    @Override
    public Class<? extends Outputer> getOutputerClass() {
        return null;
    }
}
