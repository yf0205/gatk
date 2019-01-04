package org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.segment;

public class SegmentExonOverlaps {
    private String segmentStartExonOverlap;
    private String segmentEndExonOverlap;

    public SegmentExonOverlaps(final String segmentStartExonOverlap, final String segmentEndExonOverlap) {
        this.segmentStartExonOverlap = segmentStartExonOverlap;
        this.segmentEndExonOverlap = segmentEndExonOverlap;
    }

    public String getSegmentStartExonOverlap() {
        return segmentStartExonOverlap;
    }

    public String getSegmentEndExonOverlap() {
        return segmentEndExonOverlap;
    }
}