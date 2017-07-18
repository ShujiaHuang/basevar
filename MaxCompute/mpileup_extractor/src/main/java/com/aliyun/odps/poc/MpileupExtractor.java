package com.aliyun.odps.poc;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.IllegalFormatException;
import java.util.zip.GZIPInputStream;

import com.aliyun.odps.Column;
import com.aliyun.odps.OdpsType;
import com.aliyun.odps.data.ArrayRecord;
import com.aliyun.odps.data.Record;
import com.aliyun.odps.io.InputStreamSet;
import com.aliyun.odps.io.SourceInputStream;
import com.aliyun.odps.udf.DataAttributes;
import com.aliyun.odps.udf.ExecutionContext;
import com.aliyun.odps.udf.Extractor;

/**
 * Created by lyman on 17-7-7.
 */
public class MpileupExtractor extends Extractor {

  private InputStreamSet inputs;
  private DataAttributes attributes;
  private ArrayList<String> names;
  private ArrayRecord record;
  private String line;
  private int lineNum = 0;
  private boolean firstRead = true;
  private BufferedReader mbr = null;    // buffered reader for mpileup
  private String[] tokens;
  private int nameIdx = 0;

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
  }

  private boolean moveToNextMpileup() throws IOException {
    names.clear();
    line = null;
    lineNum = 0;
    if (mbr != null) {
      mbr.close();
    }
    SourceInputStream nis = inputs.next();
    if (nis == null) {
      return false;
    }
    // read sample name file
    System.out.println("reading sample name file: " + nis.getFileName());
    BufferedReader br = new BufferedReader(new InputStreamReader(nis));
    while ((line = br.readLine()) != null) {
      names.add(line);
    }
    br.close();
    System.out.println(names.size() + " names found.");
    // open corresponding mpileup file
    SourceInputStream mis = inputs.next();
    if (mis == null) {
      System.out.println("ERR: no corresponding mpileup file!");
      return false;
    }
    System.out.println("reading mpileup file: " + mis.getFileName());
    mbr = new BufferedReader(new InputStreamReader(new GZIPInputStream(new AvailableInputStream(mis))));
    return true;
  }

  private boolean moveToNextLine() throws IOException {
    record.clear();
    nameIdx = 0;
    while (true) {
      line = mbr.readLine();
      if (line == null) {
        if (!moveToNextMpileup()) {
          return false;
        }
      } else {
        break;
      }
    }
    tokens = line.split("\t", -1);
    record.setString("chrid", tokens[0]);
    record.setString("pos", tokens[1]);
    record.setString("base_ref", tokens[2]);
    if (tokens.length != (names.size() + 1) * 3) {
      System.out.println("chrid " + tokens[0] + "\tpos " + tokens[1] + "\tbase_ref " + tokens[2]
                         + "\ttokens " + tokens.length);
      throw new IOException("mpileup data format error");
    }
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

    if (line == null || nameIdx == names.size()) {
      if (!moveToNextLine()) {
        return null;
      }
    }

    record.setString("sample_name", names.get(nameIdx++));
    int tokenIdx = nameIdx * 3; // tricky number game here
    record.setString("c1", tokens[tokenIdx++]);
    record.setString("c2", tokens[tokenIdx++]);
    record.setString("c3", tokens[tokenIdx]);
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
