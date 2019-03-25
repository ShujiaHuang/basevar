package com.aliyun.odps.exttable.ibs;

import com.aliyun.odps.udf.Extractor;
import com.aliyun.odps.udf.OdpsStorageHandler;
import com.aliyun.odps.udf.Outputer;

/**
 * Created by tianli on 2017/8/7.
 */
public class IBSHandler2 extends OdpsStorageHandler {
  public Class<? extends Extractor> getExtractorClass() {
    return IBSExtractor2.class;
  }

  public Class<? extends Outputer> getOutputerClass() {
    return null;
  }
}
