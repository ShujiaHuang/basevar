package com.aliyun.odps.exttable.ibs;

import com.aliyun.odps.udf.Extractor;
import com.aliyun.odps.udf.OdpsStorageHandler;
import com.aliyun.odps.udf.Outputer;

/**
 * Created by tianli on 2017/8/7.
 */
public class IBSHandler extends OdpsStorageHandler {
  public Class<? extends Extractor> getExtractorClass() {
    return IBSExtractor.class;
  }

  public Class<? extends Outputer> getOutputerClass() {
    return null;
  }
}
