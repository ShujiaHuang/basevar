package com.aliyun.odps.exttable.handler;

import com.aliyun.odps.exttable.output.BasevarOutputer;
import com.aliyun.odps.udf.Extractor;
import com.aliyun.odps.udf.OdpsStorageHandler;
import com.aliyun.odps.udf.Outputer;

public class BasevarHandler extends OdpsStorageHandler {
    @Override
    public Class<? extends Extractor> getExtractorClass() {
        return null;
    }

    @Override
    public Class<? extends Outputer> getOutputerClass() {
        return BasevarOutputer.class;
    }
}
