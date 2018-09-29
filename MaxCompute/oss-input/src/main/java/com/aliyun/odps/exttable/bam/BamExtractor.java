package com.aliyun.odps.exttable.bam;

import com.aliyun.odps.Column;
import com.aliyun.odps.OdpsType;
import com.aliyun.odps.data.ArrayRecord;
import com.aliyun.odps.data.Record;
import com.aliyun.odps.exttable.mpileup.AvailableInputStream;
import com.aliyun.odps.io.InputStreamSet;
import com.aliyun.odps.io.SourceInputStream;
import com.aliyun.odps.udf.DataAttributes;
import com.aliyun.odps.udf.ExecutionContext;
import com.aliyun.odps.udf.Extractor;
import htsjdk.samtools.*;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

public class BamExtractor extends Extractor {
    public final static String SAMPLE_NAME = "sample_name";
    public final static String CHRID = "chrid";
    public final static String POS = "pos";
    public final static String BASE_REF = "base_ref";
//    public final static String LINE_IDX = "line_idx";
    public final static String READ_BASE = "read_base";
    public final static String READ_QUALITY = "read_quality";
    public final static String MAPPING_QUALITY = "mapping_quality";
    public final static String READ_POS_RANK = "read_pos_rank";
    public final static String INDEL = "indel";
    public final static String STRAND = "strand";

    public final static String FASTA_FILE = "fasta_file";
    public final static String DEBUG = "debug";
    public final static String DEDUP_SIZE = "dedup_size";
    public final static String SAMPLE_FILTER = "sample_filter";
    public final static String FILE_SELECTION_MODE = "file_selection_mode";
    public final static String FILE_SELECTION_THRESHOLD = "file_selection_threshold"; // Unit: MB
    public final static String LINE_COUNT_LIMIT = "line_count_limit";


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
//    private int readQualityThreshold = 20;

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
    private String sampleName;

    private String fileSelectionMode = "small"; // 'small' or 'big_more' or 'big_less'
    private long fileSelectionThreshold = 1536l * 1024 * 1024; // Default: 1536MB
    private int lineCountLimit = 200000;


    private Set<String> sampleFilterSet = null;

    public BamExtractor() {
        Column[] columns = new Column[10];
        columns[0] = new Column(SAMPLE_NAME, OdpsType.STRING);
        columns[1] = new Column(CHRID, OdpsType.STRING);
        columns[2] = new Column(POS, OdpsType.STRING);
        columns[3] = new Column(BASE_REF, OdpsType.STRING);
//        columns[4] = new Column(LINE_IDX, OdpsType.BIGINT);
        columns[4] = new Column(READ_BASE, OdpsType.STRING);
        columns[5] = new Column(READ_QUALITY, OdpsType.STRING);
        columns[6] = new Column(MAPPING_QUALITY, OdpsType.STRING);
        columns[7] = new Column(READ_POS_RANK, OdpsType.STRING);
        columns[8] = new Column(INDEL, OdpsType.STRING);
        columns[9] = new Column(STRAND, OdpsType.STRING);
        record = new ArrayRecord(columns);
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
        String fileSelectionMode = attributes.getValueByKey(FILE_SELECTION_MODE);
        if (fileSelectionMode != null) {
            this.fileSelectionMode = fileSelectionMode;
        }
        String fileSelectionThreshold = attributes.getValueByKey(FILE_SELECTION_THRESHOLD);
        if (fileSelectionThreshold != null) {
            this.fileSelectionThreshold = Long.valueOf(fileSelectionThreshold) * 1024 * 1024;
        }
        String lineCountLimitStr = attributes.getValueByKey(LINE_COUNT_LIMIT);
        if (lineCountLimitStr != null) {
            this.lineCountLimit = Integer.valueOf(lineCountLimitStr);
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
                    sb.append(line);
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
        try {
            if (firstRead) {
                if (!moveToNextFile()) {
                    return null;
                }
                firstRead = false;
            }

            while (true) {
                if (curCigarOpratorIdx == cigarLength) {
                    if (!moveToNextCigar()) {
                        return null;
                    }
                    continue;
                }
                //curCigarOpratorIdx++; // Ensure +1 in next iteration
                switch (curCigarOpratorName) {
                    case "M":
                        curCigarOpratorIdx++;
                        int curQual = quals.charAt(curReadIdx); // char to int ascii
                        char curRead = reads.charAt(curReadIdx);
                        int rpr = curReadIdx + 1;
                        int readPos = pos;

                        curReadIdx++;
                        pos++;

                        if (deduper.isDeplicated(readPos)) {
                            // duplicate
                            continue;
                        }
                        String indelString = getIndel(); // May move cursor
//                    if (curQual - 33 >= readQualityThreshold) {
                        // Quality OK
                        record.setString(READ_BASE, String.valueOf(curRead));
                        record.setString(STRAND, strand);
                        record.setString(INDEL, indelString);
                        record.setString(READ_POS_RANK, String.valueOf(rpr)); // 1-based
                        record.setString(READ_QUALITY, String.valueOf(curQual));
                        record.setString(MAPPING_QUALITY, String.valueOf(mapq));
                        record.setString(POS, String.valueOf(readPos));
                        record.setString(BASE_REF, String.valueOf(faReference.get(chrid).charAt(readPos - 1)));
                        deduper.put(readPos);
                        return record;
//                    } else {
//                        // Quality bad
//                        if (indelString.isEmpty()) {
//                            // No indel, just skip
//                            continue;
//                        }
//                        // record indel only
//                        record.setString(READ_BASE, "N");
//                        record.setString(STRAND, ".");
//                        record.setString(INDEL, indelString);
//                        record.setString(READ_POS_RANK, "0");
//                        record.setString(READ_QUALITY, "0");
//                        record.setString(MAPPING_QUALITY, "0");
//                        record.setString(POS, String.valueOf(readPos));
//                        record.setString(BASE_REF, String.valueOf(faReference.get(chrid).charAt(readPos - 1)));
//                        deduper.put(readPos);
//                        return record;
//                    }
                    case "I":
                        curReadIdx += cigarLength;
                        curCigarOpratorIdx += cigarLength;
                        continue;
                    case "D":
                        pos += cigarLength;
                        curCigarOpratorIdx += cigarLength;
                        continue;
                    case "S": // soft
                    case "N": // jump
                    case "P": // padding
                        curReadIdx += cigarLength;
                        pos += cigarLength;
                        curCigarOpratorIdx += cigarLength;
                        continue;
                    case "H":
                        // Do nothing, just skip
                        curCigarOpratorIdx += cigarLength;
                        continue;
                    default:
                        throw new IllegalArgumentException("Unsupported cigar operator: " + curCigarOpratorName);
                }
            }
        } catch (Exception e) {
            System.err.println(String.format("chrid: %s, reads: %s, quals: %s, strand: %s, pos: %s, mapq: %s, sampleName: %s", chrid, reads, quals, strand, pos, mapq, sampleName));
            throw new RuntimeException(e);
        }
    }

    private String getIndel() throws IOException {
        if (curCigarOpratorIdx == cigarLength // last M
            && curCigarElementIdx < cigarElements.size() - 1  // has next cigar
            ) {
            moveToNextCigar();
            switch (curCigarOpratorName) {
                case "I":
                    int readStartIdx = curReadIdx;
                    curReadIdx += cigarLength;
                    curCigarOpratorIdx += cigarLength;
                    return "+" + reads.substring(readStartIdx, readStartIdx + cigarLength);
                case "D":
                    int posStartIdx = pos - 1; // pos is 1-based
                    pos += cigarLength;
                    curCigarOpratorIdx += cigarLength;
                    return "-" + faReference.get(chrid).substring(posStartIdx, posStartIdx + cigarLength);
                default:
                    return "";
            }
        }
        return "";
    }

    private boolean moveToNextFile() throws IOException {
        if (samRecordIterator != null) {
            samRecordIterator.close();
        }
        SourceInputStream nis = null;
        boolean skip;
        do {
            if (nis != null) {
                nis.close();
            }
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
            switch (fileSelectionMode) {
                case "small":
                    if (nis.getFileSize() > fileSelectionThreshold) {
                        System.err.println(String.format("[%s] File bigger than %s, skip", new Date().toLocaleString(), fileSelectionThreshold));
                        skip = true;
                    }
                    break;
                case "big_more":
                case "big_less":
                    if (nis.getFileSize() <= fileSelectionThreshold) {
                        System.err.println(String.format("[%s] File smaller than %s, skip", new Date().toLocaleString(), fileSelectionThreshold));
                        skip = true;
                    }
                    break;
                default:
                    throw new UnsupportedOperationException(fileSelectionMode);
            }
            System.out.println(String.format("[%s] Read %s, size: %s", new Date().toLocaleString(), nis.getFileName(), nis.getFileSize()));
        } while (!nis.getFileName().endsWith(".bam") || skip);

        InputStream inputStream = nis;
        if (nis.getFileName().endsWith(".gz")) {
            inputStream = new GZIPInputStream(new AvailableInputStream(nis), 65536);
        }
        SamInputResource samInputResource = SamInputResource.of(inputStream);
        SamReader samReader = srf.open(samInputResource);
        sampleName = samReader.getFileHeader().getReadGroups().get(0).getSample();
        record.setString(SAMPLE_NAME, sampleName);
        samRecordIterator = samReader.iterator();
        curLineIdx = 0;
        return true;
    }

    private boolean moveToNextLine() throws IOException {
        Set<SAMFlag> samFlags;
        String newChrId;
        boolean skipFile = false;
        while (true) {
            while (!samRecordIterator.hasNext() || skipFile) {
                if (!moveToNextFile()) {
                    return false;
                }
                skipFile = false;
            }

            // Line properties
            samRecord = samRecordIterator.next();
            curLineIdx++;
            if (curLineIdx % 1000000 == 0) {
                System.out.println(String.format("[%s] cur line: %s", new Date().toLocaleString(), curLineIdx));
            }

            if (fileSelectionMode.equalsIgnoreCase("big_more") && curLineIdx <= lineCountLimit) {
                continue;
            }
            if (fileSelectionMode.equalsIgnoreCase("big_less") && curLineIdx > lineCountLimit) {
                skipFile = true;
                continue;
            }

            newChrId = samRecord.getReferenceName();
            if(!validChridSet.contains(newChrId)) {
                if (!moveToNextFile()) {
                    return false;
                }
                continue;
            }

            samFlags = samRecord.getSAMFlags();
            mapq = samRecord.getMappingQuality();

            if (mapq >= mapqThreshold
//                && !samFlags.contains(SAMFlag.NOT_PRIMARY_ALIGNMENT)
                ) {
                break;
            }
        }
        if (!newChrId.equalsIgnoreCase(chrid)) {
            chrid = newChrId;
            deduper = new Deduper(dedupMaxSize);
        }
        pos = samRecord.getAlignmentStart();
        strand = samFlags.contains(SAMFlag.READ_REVERSE_STRAND) ? "-" : "+";

        // record set line properties
        record.setString(CHRID, chrid);
//        record.setBigint(LINE_IDX, Long.valueOf(curLineIdx));

        // Reads
        reads = samRecord.getReadString();
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
