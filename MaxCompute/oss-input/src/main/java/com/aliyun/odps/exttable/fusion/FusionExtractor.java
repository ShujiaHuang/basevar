package com.aliyun.odps.exttable.fusion;

import com.aliyun.odps.Column;
import com.aliyun.odps.OdpsType;
import com.aliyun.odps.data.ArrayRecord;
import com.aliyun.odps.data.Record;
import com.aliyun.odps.exttable.bam.Deduper;
import com.aliyun.odps.exttable.mpileup.AvailableInputStream;
import com.aliyun.odps.io.InputStreamSet;
import com.aliyun.odps.io.SourceInputStream;
import com.aliyun.odps.udf.DataAttributes;
import com.aliyun.odps.udf.ExecutionContext;
import com.aliyun.odps.udf.Extractor;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;

public class FusionExtractor extends Extractor {
    public final static String SAMPLE_NAME = "sample_name";
    public final static String CHRID = "chrid";
    public final static String START = "start_pos";
    public final static String END = "end_pos";
    public final static String ALT = "alt";
    public final static String BASE_QUALITY = "base_quality";
    public final static String READ_FIRST_POSITION = "read_first_position";
    public final static String POS = "pos";
    public final static String BASE_REF = "base_ref";
//    public final static String LINE_IDX = "line_idx";
    public final static String READ_BASE = "read_base";
    public final static String READ_QUALITY = "read_quality";
    public final static String MAPPING_QUALITY = "mapping_quality";
    public final static String READ_POS_RANK = "read_pos_rank";
    public final static String INDEL = "indel";
    public final static String STRAND = "strand";

    public final static String NON_REF = "<NON_REF>";

    public final static String FASTA_FILE = "fasta_file";
    public final static String DEBUG = "debug";
    public final static String DEDUP_SIZE = "dedup_size";
    public final static String SAMPLE_FILTER = "sample_filter";

    public final static String[] validChrids = new String[]{
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrX"
    };

    private int mapqThreshold = 30;
    private int dedupMaxSize = 2000;

    private InputStreamSet inputs;
    private ExecutionContext context;
    private SamReaderFactory samReaderFactory;
    private boolean firstRead;
    private ArrayRecord record;

    private SAMRecordIterator samRecordIterator = null;
    private SAMRecord samRecord;
    private String chrid, reads, quals, strand;
    private int pos, mapq;
    private List<CigarElement> cigarElements = null;

    private int curLineIdx = 0;
    private int curReadIdx = 0;
    private int curCigarElementIdx = 0, curCigarOpratorIdx = 0;
    private String curCigarOpratorName;
    private int cigarLength = 0;

    private Map<String, String> faReference = null;
    private Set<String> validChridSet;
    private SamReaderFactory srf;
    private Deduper deduper;
    private boolean isNewRecord = true;
    private Map<String, Integer> startPosition;
    private int curStartPosition;

    private Set<String> sampleFilterSet = null;

    public FusionExtractor() {
        Column[] columns = new Column[9];
        columns[0] = new Column(SAMPLE_NAME, OdpsType.STRING);
        columns[1] = new Column(CHRID, OdpsType.STRING);
        columns[2] = new Column(START, OdpsType.STRING);
        columns[3] = new Column(END, OdpsType.STRING);
        columns[4] = new Column(ALT, OdpsType.STRING);
        columns[5] = new Column(MAPPING_QUALITY, OdpsType.STRING);
        columns[6] = new Column(STRAND, OdpsType.STRING);
        columns[7] = new Column(READ_FIRST_POSITION, OdpsType.STRING);
        columns[8] = new Column(BASE_QUALITY, OdpsType.STRING);

        record = new ArrayRecord(columns);
        startPosition = new HashMap<>();
        validChridSet = new HashSet<>(Arrays.asList(validChrids));
    }

    @Override
    public void setup(ExecutionContext ctx, InputStreamSet inputs, DataAttributes attributes) {
        this.inputs = inputs;
        this.context = ctx;
        this.firstRead = true;
        samReaderFactory = SamReaderFactory.make();
        samReaderFactory.validationStringency(ValidationStringency.LENIENT);

        // Read reference
        String fastaFile = attributes.getValueByKey(FASTA_FILE);
        String debug = attributes.getValueByKey(DEBUG);
        if (debug == null || !debug.equalsIgnoreCase("true")) {
            readReference(fastaFile);
        }
        String dedupSize = attributes.getValueByKey(DEDUP_SIZE);
        if (dedupSize != null) {
            this.dedupMaxSize = Integer.valueOf(dedupSize);
        }
        String sampleFilter = attributes.getValueByKey(SAMPLE_FILTER);
        if (sampleFilter != null) {
            sampleFilterSet = new HashSet<>();
            try {
                BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(ctx.readResourceFileAsStream(sampleFilter)));
                String sample;
                while ((sample = bufferedReader.readLine()) != null) {
                    sampleFilterSet.add(sample);
                }
            } catch (IOException e) {
                e.printStackTrace();
                throw new RuntimeException(e);
            }
        }
        srf = SamReaderFactory.make();
        srf.validationStringency(ValidationStringency.LENIENT);
    }

    private void readReference(String fileName) {
        if (faReference != null) {
            System.out.println(String.format("[%s] Skip readReference", new Date().toLocaleString()));
            return;
        }
        faReference = new HashMap<>();
        System.out.println(String.format("[%s] Start readReference", new Date().toLocaleString()));
        BufferedReader referenceReader = null;
        try {
            BufferedInputStream bufferedInputStream = context.readResourceFileAsStream(fileName);
            referenceReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(bufferedInputStream)));
            String line;
            String curChr = null;
            StringBuilder sb = null;
            while ((line = referenceReader.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (curChr != null && sb != null && validChridSet.contains(curChr)) {
                        faReference.put(curChr, sb.toString());
                    }
                    curChr = line.substring(1).trim().split(" ")[0];
                    sb = new StringBuilder();
                } else {
                    sb.append(line.toUpperCase());
                }
            }
            if (validChridSet.contains(curChr)) {
                faReference.put(curChr, sb.toString());
            }
            System.out.println(String.format("[%s] Finish readReference. %s",
                new Date().toLocaleString(),
                faReference.keySet()));
        } catch (IOException e) {
            throw new RuntimeException("Read fasta file failed!", e);
        } finally {
            if (referenceReader != null) {
                try {
                    referenceReader.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    @Override
    public Record extract() throws IOException {
        if (firstRead) {
            if (!moveToNextFile()) {
                return null;
            }
            firstRead = false;
        }

        while (true) {
            if (curCigarOpratorIdx == cigarLength) {
                if (!isNewRecord) {
                    isNewRecord = true;
                    return record;
                }
                if (!moveToNextCigar()) {
                    return null;
                }
                continue;
            }
            if (pos < curStartPosition) {
                curCigarOpratorIdx++;
                pos++;
                curReadIdx++;
                continue;
            }
            switch (curCigarOpratorName) {
                case "M":
                    char curQual = quals.charAt(curReadIdx); // char to int ascii
                    char curRead = reads.charAt(curReadIdx);
                    char curRefBase = faReference.get(chrid).charAt(pos);

                    if (isNewRecord) {
                        isNewRecord = false;
                        record.setString(START, String.valueOf(pos));
                        record.setString(END, String.valueOf(pos + 1));
                        record.setString(READ_FIRST_POSITION, String.valueOf(curReadIdx + 1));
                        record.setString(BASE_QUALITY, String.valueOf(curQual));
                        record.setString(ALT, (curRead == curRefBase) ? NON_REF : String.valueOf(curRead));
                        curReadIdx++;
                        pos++;
                        curCigarOpratorIdx++;
                    } else if (record.getString(ALT).equals(NON_REF)
                        && curRead == curRefBase
                        && record.getString(BASE_QUALITY).charAt(0) == curQual) {
                        record.setString(END, String.valueOf(pos + 1));
                        curCigarOpratorIdx++;
                        curReadIdx++;
                        pos++;
                    } else {
                        isNewRecord = true;
                        return record;
                    }
                    break;
                case "I":
                    if (!isNewRecord) {
                        isNewRecord = true;
                        return record;
                    }
                    record.setString(START, String.valueOf(pos - 1));
                    record.setString(END, String.valueOf(pos));
                    record.setString(ALT, reads.substring(curReadIdx, curReadIdx+cigarLength));
                    record.setString(READ_FIRST_POSITION, "0"); // No meaning
                    record.setString(BASE_QUALITY, quals.substring(curReadIdx,curReadIdx+cigarLength));
                    curReadIdx += cigarLength;
                    curCigarOpratorIdx += cigarLength;
                    return record;
                case "D":
                    if (!isNewRecord) {
                        isNewRecord = true;
                        return record;
                    }
                    record.setString(START, String.valueOf(pos));
                    record.setString(END, String.valueOf(pos+cigarLength));
                    record.setString(ALT, ".");
                    record.setString(READ_FIRST_POSITION, "0"); // No meaning
                    record.setString(BASE_QUALITY, "!");
                    pos += cigarLength;
                    curCigarOpratorIdx += cigarLength;
                    return record;
                case "S": // soft
                case "N": // jump
                case "P": // padding
                    curReadIdx += cigarLength;
                    pos += cigarLength;
                    curCigarOpratorIdx += cigarLength;
                    if (!isNewRecord) {
                        isNewRecord = true;
                        return record;
                    }
                    break;
                case "H":
                    // Do nothing, just skip
                    curCigarOpratorIdx += cigarLength;
                    if (!isNewRecord) {
                        isNewRecord = true;
                        return record;
                    }
                    break;
                default:
                    throw new IllegalArgumentException("Unsupported cigar operator: " + curCigarOpratorName);
            }
        }
    }

    private boolean moveToNextFile() throws IOException {
        if (samRecordIterator != null) {
            samRecordIterator.close();
        }
        SourceInputStream nis;
        boolean skip;
        do {
            nis = inputs.next();
            if (nis == null) {
                return false;
            }
            skip = false;
            if (sampleFilterSet != null) {
                String[] tokens = nis.getFileName().split("/");
                if (!sampleFilterSet.contains(tokens[tokens.length - 1].split("\\.")[0])) {
                    System.out.println(String.format("[%s] Skip %s", new Date().toLocaleString(), nis.getFileName()));
                    skip = true;
                    continue;
                }
            }
            System.out.println(String.format("[%s] Read %s", new Date().toLocaleString(), nis.getFileName()));
        } while (!nis.getFileName().endsWith(".bam") || skip);

        InputStream inputStream = new AvailableInputStream(nis);
        if (nis.getFileName().endsWith(".gz")) {
            inputStream = new GZIPInputStream(inputStream, 65536);
        }
        SamInputResource samInputResource = SamInputResource.of(inputStream);
        SamReader samReader = srf.open(samInputResource);
        String sampleName = samReader.getFileHeader().getReadGroups().get(0).getSample();
        record.setString(SAMPLE_NAME, sampleName);
        samRecordIterator = samReader.iterator();
        curLineIdx = 0;
        return true;
    }

    private boolean moveToNextLine() throws IOException {
        Set<SAMFlag> samFlags;
        String newChrId;
        if (samRecord != null &&
            samRecord.getAlignmentEnd() > startPosition.getOrDefault(samRecord.getReferenceName(), 0)) {
            startPosition.put(samRecord.getReferenceName(), samRecord.getAlignmentEnd());
        }
        while (true) {
            while (!samRecordIterator.hasNext()) {
                if (!moveToNextFile()) {
                    return false;
                }
            }
            // Line properties
            samRecord = samRecordIterator.next();
            newChrId = samRecord.getReferenceName();
            if(!validChridSet.contains(newChrId)) {
                if (!moveToNextFile()) {
                    return false;
                }
                continue;
            }
            curLineIdx++;
            if (curLineIdx % 1000000 == 0) {
                System.out.println(String.format("[%s] cur line: %s", new Date().toLocaleString(), curLineIdx));
            }
            samFlags = samRecord.getSAMFlags();
            mapq = samRecord.getMappingQuality();

            if (mapq >= mapqThreshold) {
                break;
            }
        }
        curStartPosition = startPosition.getOrDefault(samRecord.getReferenceName(), 0);
        if (!newChrId.equalsIgnoreCase(chrid)) {
            chrid = newChrId;
            deduper = new Deduper(dedupMaxSize);
        }
        pos = samRecord.getAlignmentStart() - 1;
        strand = samFlags.contains(SAMFlag.READ_REVERSE_STRAND) ? "-" : "+";

        // record set line properties
        record.setString(CHRID, chrid);
        record.setString(MAPPING_QUALITY, String.valueOf(mapq));
        record.setString(STRAND, strand);
//        record.setBigint(LINE_IDX, Long.valueOf(curLineIdx));

        // Reads
        reads = samRecord.getReadString().toUpperCase();
        quals = samRecord.getBaseQualityString();
        cigarElements = samRecord.getCigar().getCigarElements();

        // Idx
        curCigarElementIdx = 0; // ++ in move to next cigar
        curReadIdx = 0;

        return true;
    }

    private boolean moveToNextCigar() throws IOException {
        if (cigarElements == null || curCigarElementIdx == cigarElements.size()) {
            if (!moveToNextLine()) {
                return false;
            }
        }

        curCigarOpratorIdx = 0;
        CigarElement cigarElement = cigarElements.get(curCigarElementIdx);
        curCigarOpratorName = cigarElement.getOperator().name();
        cigarLength = cigarElement.getLength();
        curCigarElementIdx++;
        return true;
    }

    @Override
    public void close() {
        if (samRecordIterator != null) {
            samRecordIterator.close();
        }
    }
}
