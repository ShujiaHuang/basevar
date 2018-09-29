package com.aliyun.odps.exttable.mpileup;

import com.aliyun.odps.io.SourceInputStream;

import java.io.IOException;
import java.nio.BufferOverflowException;

/**
 * This class is a workaround for http://bugs.java.com/view_bug.do?bug_id=7021870
 * by always returning 1 at available()
 * Created by lyman on 17-7-18.
 */
public class AvailableInputStream extends SourceInputStream {

  private SourceInputStream sis;

  public AvailableInputStream(SourceInputStream is) {
    this.sis = is;
  }

  public int read(byte[] bytes, int i, int i1) throws IOException {
    return sis.read(bytes, i, i1);
  }

  public int read(byte[] bytes) throws IOException {
    return sis.read(bytes);
  }

  public int read() throws IOException {
    return sis.read();
  }

  public String getFileName() {
    return sis.getFileName();
  }

  public long getFileSize() {
    return sis.getFileSize();
  }

  public int readToEnd(byte[] bytes) throws IOException, BufferOverflowException {
    return sis.readToEnd(bytes);
  }

  @Override
  public int available() throws IOException {
    return 1; // always return available
  }
}
