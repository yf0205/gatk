package org.broadinstitute.hellbender.utils.downsampling;

import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.commons.lang.mutable.MutableInt;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;


/**
 * PositionalDownsampler: Downsample each stack of reads at each alignment start to a size <= a target coverage
 * using a {@link ReservoirDownsampler}. Stores only O(target coverage) reads in memory at any given time,
 * provided the client regularly calls {@link #consumeFinalizedItems}.
 *
 * Unmapped reads with assigned positions are subject to downsampling in the same way as mapped reads,
 * but unmapped reads without assigned positions are not subject to downsampling.
 */
public final class MutectDownsampler extends ReadsDownsampler {

    // we reject reads with a lot of low-base quality substitutions.  These reads can't be trusted.
    private static final int BAD_BASE_QUALITY = 15;
    private final int maxLowQualSubstitutionsInRead;

    // we don't reject individual reads based on a lot of substitutions or mediocre mapping quality,
    // but if a lot of reads are suspicious in this way (i.e. abnormally high coverage) we reject the whole locus
    // as a mapping error hotspot
    private static final int SUSPICIOUS_NUM_MISMATCHES = 3;
    private static final int SUSPICIOUS_MAPPING_QUALITY = 55;
    private final MutableInt suspiciousReadCount;
    private final int maxSuspiciousReadsPerStride;

    private final int stride;
    private final int maxCoverage;
    private final List<GATKRead> pendingReads;

    private final SAMFileHeader header;
    private GATKRead firstReadInStride;
    private List<GATKRead> finalizedReads;
    private final ReferenceDataSource referenceDataSource;

    private byte[] refBases;
    private int refBasesStart = -1;

    private boolean rejectAllReadsInStride;

    // when we fetch ref bases, add a buffer in order to reuse these bases for more reads
    private static final int REF_BASES_EXTRA_BUFFER = 10000;

    /**
     * Construct a PositionalDownsampler
     *
     * @param maxReadsPerAlignmentStart Maximum number of reads per alignment start position. Must be > 0
     * @param stride Length in bases constituting a single pool of reads to downsample
     * @param header SAMFileHeader to use to determine contig ordering. Non-null.
     */
    public MutectDownsampler(final int maxReadsPerAlignmentStart,
                             final SAMFileHeader header,
                             final int maxLowQualSubstitutionsInRead,
                             final int maxSuspiciousReadsPerAlignmentStart,
                             final int stride,
                             final ReferenceDataSource referenceDataSource) {
        this.header = Utils.nonNull(header);

        // convert coverage per base to coverage per stride
        maxCoverage = maxReadsPerAlignmentStart <= 0 ? Integer.MAX_VALUE : maxReadsPerAlignmentStart * stride;
        this.stride = ParamUtils.isPositive(stride, "stride must be > 0");
        this.maxLowQualSubstitutionsInRead = ParamUtils.isPositive(maxLowQualSubstitutionsInRead, "maxLowQualSubstitutionsInRead must be > 0");
        maxSuspiciousReadsPerStride = stride * ParamUtils.isPositive(maxSuspiciousReadsPerAlignmentStart, "maxSuspiciousReadsPerAlignmentStart must be > 0");
        this.referenceDataSource = referenceDataSource;

        pendingReads = new ArrayList<>();
        finalizedReads = new ArrayList<>();
        rejectAllReadsInStride = false;
        suspiciousReadCount = new MutableInt(0);

        clearItems();
        resetStats();
    }

    @Override
    public void submit( final GATKRead newRead ) {
        Utils.nonNull(newRead);
        if (ReadUtils.readHasNoAssignedPosition(newRead)) {
            finalizedReads.add(newRead);
            return;
        }

        handlePositionalChange(newRead);

        if (rejectAllReadsInStride) {
            return;
        }

        final Pair<Integer, Integer> highQualAndLowQualMismatches = countMismatches(newRead);
        final int highQualMismatches = highQualAndLowQualMismatches.getLeft();
        final int lowQualMismatches = highQualAndLowQualMismatches.getRight();

        if (highQualMismatches >= SUSPICIOUS_NUM_MISMATCHES || newRead.getMappingQuality() <= SUSPICIOUS_MAPPING_QUALITY) {
            suspiciousReadCount.increment();
        }

        rejectAllReadsInStride |= suspiciousReadCount.intValue() >= maxSuspiciousReadsPerStride;

        if (lowQualMismatches <= maxLowQualSubstitutionsInRead) {
            pendingReads.add(newRead);
        }
    }

    // count the number of non-consecutive single-base mismatches between a read and the reference
    // soft clips and indels are not counted, and MNPs are counted as one mismatch.
    // high and low base quality mismatches are counted separately because the former indicates possible
    // mapping error, while the latter indicates sequencer error.
    private Pair<Integer, Integer> countMismatches(final GATKRead read) {
        // see if we can use an existing NM SAM tag
        if (read.getCigarElements().size() == 1) {
            final Integer editDistance = read.getAttributeAsInteger(SAMTag.NM.name());
            if (editDistance == 0 || editDistance == 1) {
                return new ImmutablePair<>(editDistance, 0);
            }
        }

        // converting to SAM record avoids defensive copying when getting bases and quals
        final SAMRecord read1 = read.convertToSAMRecord(header);
        final MutableInt highQualMismatches = new MutableInt(0);
        final MutableInt lowQualMismatches = new MutableInt(0);
        final byte[] readBases = read1.getReadBases();
        final byte[] readBaseQualities = read1.getBaseQualities();

        for (final AlignmentBlock block : read1.getAlignmentBlocks()) {
            final int readBlockStart = block.getReadStart() - 1;
            final int referenceBlockStart = block.getReferenceStart() - refBasesStart;
            final int length = block.getLength();

            boolean lastBaseWasMatch = true;
            for (int i = 0; i < length; ++i) {
                final byte readBase = readBases[readBlockStart + i];
                final int referenceIndex = referenceBlockStart + i;
                final byte refBase = refBases[referenceIndex];
                final boolean basesMatch = SequenceUtil.basesEqual(readBase, refBase);
                if (!basesMatch && lastBaseWasMatch) {
                    final byte readBaseQuality = readBaseQualities[readBlockStart + i];
                    if (readBaseQuality > BAD_BASE_QUALITY) {
                        highQualMismatches.increment();
                    } else {
                        lowQualMismatches.increment();
                    }
                }
                lastBaseWasMatch = basesMatch;
            }
        }

        return new ImmutablePair<>(highQualMismatches.intValue(), lowQualMismatches.intValue());

    }

    private void handlePositionalChange( final GATKRead newRead ) {
        if (ReadUtils.readHasNoAssignedPosition(newRead)) {
            return;
        } else if (firstReadInStride == null) {
            firstReadInStride = newRead;
            fetchRefBases(newRead);
        } else {
            final boolean newContig = !newRead.getAssignedContig().equals(firstReadInStride.getAssignedContig());
            final boolean newStride = newContig || newRead.getAssignedStart() >= firstReadInStride.getAssignedStart() + stride;
            final boolean needsNewRefBases = newContig || newRead.getEnd() >= refBasesStart + refBases.length - 1;
            if (newStride) {
                finalizePendingReads();
                firstReadInStride = newRead;
                rejectAllReadsInStride = false;
                suspiciousReadCount.setValue(0);
            }

            if (needsNewRefBases) {
                fetchRefBases(newRead);
            }
        }
    }

    private void fetchRefBases(GATKRead newRead) {
        refBases = referenceDataSource.queryAndPrefetch(newRead.getAssignedContig(), newRead.getAssignedStart(), newRead.getEnd() + REF_BASES_EXTRA_BUFFER).getBases();
        refBasesStart = newRead.getAssignedStart();
    }


    private void finalizePendingReads() {
        if (!rejectAllReadsInStride) {
            // most common case: no downsampling necessary
            if (pendingReads.size() <= maxCoverage) {
                finalizedReads.addAll(pendingReads);
            } else {
                // if we exceed the max coverage, just use well-mapped reads.  Maybe the number of such reads won't reach
                // the desired coverage, but if the region is decently mappable the shortfall will be minor.
                final ReservoirDownsampler wellMappedDownsampler = new ReservoirDownsampler(maxCoverage, false);
                pendingReads.stream().filter(read -> read.getMappingQuality() > SUSPICIOUS_MAPPING_QUALITY).forEach(wellMappedDownsampler::submit);
                final List<GATKRead> readsToFinalize = wellMappedDownsampler.consumeFinalizedItems();
                if (stride > 1) {
                    Collections.sort(readsToFinalize, Comparator.comparingInt(GATKRead::getAssignedStart));
                }
                finalizedReads.addAll(readsToFinalize);
            }
        }
        pendingReads.clear();

    }

    @Override
    public boolean hasFinalizedItems() {
        return ! finalizedReads.isEmpty();
    }

    @Override
    public boolean hasPendingItems() { return !pendingReads.isEmpty(); }

    @Override
    public GATKRead peekFinalized() { return finalizedReads.isEmpty() ? null : finalizedReads.get(0); }

    @Override
    public GATKRead peekPending() { return pendingReads.isEmpty() ? null : pendingReads.get(0); }

    @Override
    public List<GATKRead> consumeFinalizedItems() {
        // if stride is > 1, different starts may have gotten mixed up
        // TODO: this might not end up being relevant
        if (stride > 1) {
            Collections.sort(finalizedReads, Comparator.comparingInt(GATKRead::getAssignedStart));
        }
        final List<GATKRead> toReturn = finalizedReads;
        finalizedReads = new ArrayList<>();
        return toReturn;
    }

    @Override
    public int size() { return finalizedReads.size() + pendingReads.size(); }

    @Override
    public void signalEndOfInput() {
        finalizePendingReads();
    }

    @Override
    public void clearItems() {
        pendingReads.clear();
        finalizedReads.clear();
        firstReadInStride = null;
        rejectAllReadsInStride = false;
        suspiciousReadCount.setValue(0);
    }

    @Override
    public boolean requiresCoordinateSortOrder() {
        return true;
    }

    @Override
    public void signalNoMoreReadsBefore( final GATKRead read ) {
        handlePositionalChange(read);
    }
}
