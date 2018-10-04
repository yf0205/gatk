package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.List;
import java.util.stream.Collectors;

public class AlleleFractionClustering {
    private AlleleFractionClusterer somaticClustering;
    private double log10SomaticPrior;
    private double log10NothingPrior;
    private double posteriorThreshold;

    private static final double CONCENTRATION = 0.1;

    public AlleleFractionClustering(final List<ImmutablePair<double[], int[]>> tumorLodsAndCounts,
                                    final long callableSites, final M2FiltersArgumentCollection MTFAC) {

        final List<ImmutablePair<Double, Count>> lodsAndCounts = tumorLodsAndCounts.stream().map(pair -> {
            final double bestTLod = MathUtils.arrayMax(pair.getLeft());
            final int refCount = pair.getRight()[0];
            final int altCount = pair.getRight()[1];    //TODO: handle multiallelics
            return ImmutablePair.of(bestTLod, new Count(altCount, refCount));
        }).collect(Collectors.toList());

        final List<Count> confidentCounts = lodsAndCounts.stream()
                .filter(pair -> pair.getLeft() > MTFAC.highConfidenceLod)
                .map(pair -> pair.getRight())
                .collect(Collectors.toList());

        somaticClustering = new AlleleFractionClusterer(confidentCounts, CONCENTRATION);

        final double somaticCount = confidentCounts.size();
        log10SomaticPrior = FastMath.log10(somaticCount / callableSites);
        log10NothingPrior = FastMath.log10((callableSites - somaticCount)/ callableSites);

        // TODO: put this in convergence loop where we update the priors?
        final double[] posteriorProbabilitiesOfNothing = lodsAndCounts.stream()
                .mapToDouble(pair -> 1 - somaticProbability(pair.getLeft(), pair.getRight())).toArray();

        final FilteringFirstPass.FilterStats filterStats = FilteringFirstPass.calculateFilterThreshold(posteriorProbabilitiesOfNothing, MTFAC.maxFalsePositiveRate, GATKVCFConstants.TUMOR_LOD_FILTER_NAME);
        posteriorThreshold = filterStats.getThreshold();
    }

    // array of responsibilities / posteriors, starting with nothing and proceeding through all somatic clusters
    // (currently there is only one somatic cluster)
    private double somaticProbability(final double tumorLog10Odds, final Count count) {
        final double[] unweightedLog10Responsibilities = new double[]{log10NothingPrior,
                log10SomaticPrior + tumorLog10Odds + somaticClustering.log10OddsCorrection(count.getAltCount(), count.getRefCount())};

        return MathUtils.normalizeFromLog10ToLinearSpace(unweightedLog10Responsibilities)[1];
    }

    public boolean passesThreshold(final double tumorLog10Odds, final int refCount, final int altCount) {
        final double nonSomaticProbability = 1 - somaticProbability(tumorLog10Odds, new Count(altCount, refCount));
        return nonSomaticProbability < posteriorThreshold;
    }

}
