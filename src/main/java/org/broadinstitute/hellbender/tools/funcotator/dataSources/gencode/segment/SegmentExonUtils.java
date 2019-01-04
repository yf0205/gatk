package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.segment;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.annotation.Strand;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfExonFeature;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfTranscriptFeature;

import java.util.List;

public class SegmentExonUtils {
    private SegmentExonUtils(){}

    // This should always be -1
    private final static int NO_EXON_OVERLAP = -1;

    public static SegmentExonOverlaps determineSegmentExonPosition(final GencodeGtfTranscriptFeature transcript, final Locatable segment, final SAMSequenceDictionary dictionary) {
        final List<GencodeGtfExonFeature> exons = transcript.getExons();
        final String NOT_IN_TRANSCRIPT = "";

        if (exons.size() == 0) {
            return new SegmentExonOverlaps(NOT_IN_TRANSCRIPT, NOT_IN_TRANSCRIPT);
        }

        final boolean[] isExonOverlappedMask = createExonOverlapMask(segment, exons);
        final SimpleInterval segmentStart = new SimpleInterval(segment.getContig(), segment.getStart(), segment.getStart());
        final SimpleInterval segmentEnd = new SimpleInterval(segment.getContig(), segment.getEnd(), segment.getEnd());

        // Find the proper index for the start.
        // If the start of the segment does not overlap the first exon, but the segment does, then the start index is -1 (no overlap)
        final int inclusiveIndexPositiveDirectionStart = findInclusiveIndexPositiveDirectionStart(transcript, exons, isExonOverlappedMask, segmentStart);

        // If the end of the segment does not overlap the last exon, but the segment does, then the end index is -1 (no overlap)
        final int inclusiveIndexPositiveDirectionEnd = findInclusiveIndexPositiveDirectionEnd(transcript, exons, isExonOverlappedMask, segmentEnd);

        // Construct the final strings
        final String startResult = inclusiveIndexPositiveDirectionStart != NO_EXON_OVERLAP ?
                inclusiveIndexPositiveDirectionStart + determineSegmentOverlapDirection(transcript.getGenomicStrand(), true)
                : NOT_IN_TRANSCRIPT;
        final String endResult = inclusiveIndexPositiveDirectionEnd != NO_EXON_OVERLAP ?
                inclusiveIndexPositiveDirectionEnd + determineSegmentOverlapDirection(transcript.getGenomicStrand(), false)
                : NOT_IN_TRANSCRIPT;

        return new SegmentExonOverlaps(startResult, endResult);
    }

    private static int findInclusiveIndexPositiveDirectionEnd(final GencodeGtfTranscriptFeature transcript, final List<GencodeGtfExonFeature> exons, final boolean[] isExonOverlappedMask, final SimpleInterval segmentEnd) {
        int inclusiveIndexPositiveDirectionEnd = NO_EXON_OVERLAP;
        final int lastExonIndex = exons.size() - 1;
        if (isExonOverlappedMask[lastExonIndex] && IntervalUtils.overlaps(segmentEnd, exons.get(lastExonIndex))) {
            inclusiveIndexPositiveDirectionEnd = convertInclusiveIndexForNegativeStrand(transcript, lastExonIndex);
        } else if (!isExonOverlappedMask[lastExonIndex]) {

            // Find the first exon that has overlap with the segment, going backwards in genomic space
            for (int i = lastExonIndex - 1; i >= 0; i --) {
                if (isExonOverlappedMask[i]) {
                    inclusiveIndexPositiveDirectionEnd =  convertInclusiveIndexForNegativeStrand(transcript, i);
                    break;
                }
            }
        } else {
            inclusiveIndexPositiveDirectionEnd = NO_EXON_OVERLAP;
        }
        return inclusiveIndexPositiveDirectionEnd;
    }

    private static int findInclusiveIndexPositiveDirectionStart(final GencodeGtfTranscriptFeature transcript, final List<GencodeGtfExonFeature> exons, final boolean[] isExonOverlappedMask, final SimpleInterval segmentStart) {
        int inclusiveIndexPositiveDirectionStart = NO_EXON_OVERLAP;
        if (isExonOverlappedMask[0] && IntervalUtils.overlaps(segmentStart, exons.get(0))) {
            inclusiveIndexPositiveDirectionStart = convertInclusiveIndexForNegativeStrand(transcript,0);
        } else if (!isExonOverlappedMask[0]) {
            // Find the first exon that has overlap with the segment
            for (int i = 1; i < exons.size(); i ++) {
                if (isExonOverlappedMask[i]) {
                    inclusiveIndexPositiveDirectionStart =  convertInclusiveIndexForNegativeStrand(transcript, i);
                    break;
                }
            }
        } else {
            inclusiveIndexPositiveDirectionStart = NO_EXON_OVERLAP;
        }
        return inclusiveIndexPositiveDirectionStart;
    }

    private static boolean[] createExonOverlapMask(final Locatable segment, final List<GencodeGtfExonFeature> exons) {
        final boolean[] isExonOverlappedMask = new boolean[exons.size()];

        for (int i = 0; i < exons.size(); i++) {
            final GencodeGtfExonFeature exon = exons.get(i);

            if (IntervalUtils.overlaps(segment, exon.getGenomicPosition())) {
                isExonOverlappedMask[i] = true;
            }
        }
        return isExonOverlappedMask;
    }

    @VisibleForTesting
    static String determineSegmentOverlapDirection(final Strand strand, final boolean isSegmentStart) {
        if (isSegmentStart ^ (strand == Strand.POSITIVE)) {
            return "-";
        } else {
            return "+";
        }
    }

    @VisibleForTesting
    static int convertInclusiveIndexForNegativeStrand(final GencodeGtfTranscriptFeature transcript, final int initialInclusiveIndexPositiveDirection) {
        final List<GencodeGtfExonFeature> exons = transcript.getExons();

        if ((transcript.getGenomicStrand() == Strand.NEGATIVE) && (initialInclusiveIndexPositiveDirection != NO_EXON_OVERLAP)){
            return exons.size() - initialInclusiveIndexPositiveDirection - 1;
        }
        return initialInclusiveIndexPositiveDirection;
    }
}
