package org.broadinstitute.hellbender.tools.walkers.mutect;

import java.util.List;
import java.util.Set;
import java.util.stream.IntStream;

public class AlleleFractionClusterer {

    final int INITIAL_CRP_CLUSTERS = 20;

    public AlleleFractionClusterer(final List<Count> altAndRefCounts) {
        final double[] sortedFractions = altAndRefCounts.stream().mapToDouble(p -> p.getAF()).sorted().toArray();

        final int numInitialClusters = Math.min(altAndRefCounts.size(), INITIAL_CRP_CLUSTERS);

        final double[] initialClusters = IntStream.range(1, numInitialClusters + 1).mapToDouble(n -> {
            final double percentile = ((double) n / numInitialClusters) * 100;
            return
        })

    }

    public static class Count {
        int altCount;
        int refCount;

        public Count(final int altCount, final int refCount) {
            this.altCount = altCount;
            this.refCount = refCount;
        }

        public double getAF() { return (double) altCount / refCount; }
    }

    private static class Cluster {
        double alleleFraction;
        final Set<Count> counts;
    }
}
