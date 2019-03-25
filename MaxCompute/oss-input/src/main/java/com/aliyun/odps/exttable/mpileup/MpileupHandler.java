package com.aliyun.odps.exttable.mpileup;

import com.aliyun.odps.udf.Extractor;
import com.aliyun.odps.udf.OdpsStorageHandler;
import com.aliyun.odps.udf.Outputer;

/**
 * Created by lyman on 17-7-7.
 */
public class MpileupHandler extends OdpsStorageHandler {

  @Override
  public Class<? extends Extractor> getExtractorClass() {
    return MpileupExtractor.class;
  }
  @Override
  public Class<? extends Outputer> getOutputerClass() {
    return null;
  }

}
