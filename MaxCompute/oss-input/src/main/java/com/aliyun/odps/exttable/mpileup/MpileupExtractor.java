package com.aliyun.odps.exttable.mpileup;

import com.aliyun.odps.Column;
import com.aliyun.odps.OdpsType;
import com.aliyun.odps.data.ArrayRecord;
import com.aliyun.odps.data.Record;
import com.aliyun.odps.io.InputStreamSet;
import com.aliyun.odps.io.SourceInputStream;
import com.aliyun.odps.udf.DataAttributes;
import com.aliyun.odps.udf.ExecutionContext;
import com.aliyun.odps.udf.Extractor;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.zip.GZIPInputStream;

/**
 * This process one-sample mpile files
 */
public class MpileupExtractor extends Extractor {

  private InputStreamSet inputs;
  private DataAttributes attributes;
  private ArrayList<String> names;
  private ArrayRecord record;
  private String line;
  private boolean firstRead = true;
  private BufferedReader mbr = null;    // buffered reader for mpileup
  private String[] tokens;
  private String curSampleName = "";
  private ExecutionContext context;

  public MpileupExtractor() {
    Column[] columns = new Column[7];
    columns[0] = new Column("chrid", OdpsType.STRING);
    columns[1] = new Column("pos", OdpsType.STRING);
    columns[2] = new Column("base_ref", OdpsType.STRING);
    columns[3] = new Column("sample_name", OdpsType.STRING);
    columns[4] = new Column("c1", OdpsType.STRING);
    columns[5] = new Column("c2", OdpsType.STRING);
    columns[6] = new Column("c3", OdpsType.STRING);
    record = new ArrayRecord(columns);
    names = new ArrayList<String>();
  }

  // no particular usage for execution context in this example
  @Override
  public void setup(ExecutionContext ctx, InputStreamSet inputs, DataAttributes attributes) {
    this.inputs = inputs;
    this.attributes = attributes;
    this.context = ctx;
  }

  private boolean moveToNextMpileup() throws IOException {
    names.clear();
    line = null;
    if (mbr != null) {
      mbr.close();
    }
    SourceInputStream nis = inputs.next();
    if (nis == null) {
      return false;
    }
    // read sample name file
    String fileName = nis.getFileName();
    String[] fileArray = fileName.split("\\.")[0].split("/");
    curSampleName = fileArray[fileArray.length - 1];
    record.setString("sample_name", curSampleName);
    System.out.println("reading sample name file: " + fileName);

    mbr = new BufferedReader(new InputStreamReader(new GZIPInputStream(new AvailableInputStream(nis))));
    return true;
  }

  private boolean moveToNextLine() throws IOException {
    while (true) {
      line = mbr.readLine();

      if (line == null) {
        if (!moveToNextMpileup()) {
          return false;
        }
      } else {
        tokens = line.split("\t", -1);
        if (tokens[3].equalsIgnoreCase("0")
            && tokens[4].equalsIgnoreCase("*")
            && tokens[5].equalsIgnoreCase("*")) {
          continue;
        }
        break;
      }
    }

    if (tokens.length != 6) {
      System.out.println("chrid " + tokens[0] + "\tpos " + tokens[1] + "\tbase_ref " + tokens[2]
          + "\ttokens " + tokens.length);
      throw new IOException("mpileup data format error");
    }

    record.setString("chrid", tokens[0]);
    record.setString("pos", tokens[1]);
    record.setString("base_ref", tokens[2]);
    record.setString("c1", tokens[3]);
    record.setString("c2", tokens[4]);
    record.setString("c3", tokens[5]);

    return true;
  }

  @Override
  public Record extract() throws IOException {
    if (firstRead) {
      if (!moveToNextMpileup()) {
        return null;
      }
      firstRead = false;
    }

    if (!moveToNextLine()) {
      return null;
    }
    return record;
  }

  @Override
  public void close(){
    // no-op
  }

  public static void main(String[] args) throws IOException {

    FileInputStream fis = new FileInputStream(args[0]);
    BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(fis)));
    String line;
    while ((line = br.readLine()) != null) {
      String[] tokens = line.split("\t", -1);
      System.out.println(tokens[0] + "\t" + tokens[1] + "\t" + tokens[2] + "\t" + tokens.length);
    }
  }
}
