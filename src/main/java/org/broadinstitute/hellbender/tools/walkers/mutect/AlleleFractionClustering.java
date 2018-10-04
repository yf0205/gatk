package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.List;
import java.util.stream.Collectors;

public class AlleleFractionClustering {
    private AlleleFractionClusterer somaticDistribution;
    private final long callableSites;
    private double log10SomaticPrior;
    private double log10NothingPrior;
    private double posteriorThreshold;

    private static final double CONCENTRATION = 0.1;

    public AlleleFractionClustering(final List<ImmutablePair<double[], int[]>> tumorLodsAndCounts,
                                    final long callableSites, final M2FiltersArgumentCollection MTFAC) {
        this.callableSites = callableSites;

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

        somaticDistribution = new AlleleFractionClusterer(confidentCounts, CONCENTRATION);

        final double somaticCount = confidentCounts.size();
        log10SomaticPrior = FastMath.log10(somaticCount / callableSites);
        log10NothingPrior = FastMath.log10((callableSites - somaticCount)/ callableSites);

        //TODO: real measure of convergence or at least extract magic constant!
        //for (int n = 0; n < 3; n++) {
            // E step with fractional assignments to nothing, low-confidence, high confidence
            responsibilities = tumorLodsAndCounts.stream().map(pair -> {
                final double tumorLog10Odds = pair.getLeft()[0];
                final double refCount = pair.getRight()[0];
                final double altCount = pair.getRight()[1];
                return getResponsibilities(tumorLog10Odds, refCount, altCount);
            }).collect(Collectors.toList());

        //    updatePriors(responsibilities);
        //    fitShape(allCounts, responsibilities);
        //}

        final double[] posteriorProbabilitiesOfNothing = responsibilities.stream().mapToDouble(r -> r[0]).toArray();
        final FilteringFirstPass.FilterStats filterStats = FilteringFirstPass.calculateFilterThreshold(posteriorProbabilitiesOfNothing, MTFAC.maxFalsePositiveRate, GATKVCFConstants.TUMOR_LOD_FILTER_NAME);
        posteriorThreshold = filterStats.getThreshold();
    }

    // array of responsibilities / posteriors, starting with nothing and proceeding through all somatic clusters
    // (currently there is only one somatic cluster)
    private double[] getResponsibilities(final double tumorLog10Odds, final double refCount, final double altCount) {
        final double log10OddsCorrection = SomaticLikelihoodsEngine.log10OddsCorrection(
                somaticDistribution, Dirichlet.flat(2), new double[]{altCount, refCount});

        final double[] unweightedLog10Responsibilities = new double[]{log10NothingPrior,
                log10SomaticPrior + tumorLog10Odds + log10OddsCorrection};

        return MathUtils.normalizeFromLog10ToLinearSpace(unweightedLog10Responsibilities);
    }

    public boolean passesThreshold(final double tumorLog10Odds, final double refCount, final double altCount) {
        final double nonSomaticProbability =  getResponsibilities(tumorLog10Odds, refCount, altCount)[0];
        return nonSomaticProbability < posteriorThreshold;
    }

}
