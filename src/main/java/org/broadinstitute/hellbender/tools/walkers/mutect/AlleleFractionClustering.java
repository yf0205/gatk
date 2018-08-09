package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math.MathException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.solvers.BisectionSolver;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class AlleleFractionClustering {
    private BetaDistributionShape somaticDistribution = new BetaDistributionShape(1,1);
    private final long callableSites;
    private double log10SomaticPrior;
    private double log10NothingPrior;
    private double posteriorThreshold;

    public AlleleFractionClustering(final List<ImmutablePair<double[], double[]>> tumorLodsAndCounts,
                                    final long callableSites, final M2FiltersArgumentCollection MTFAC) {
        this.callableSites = callableSites;
        final List<double[]> allCounts = tumorLodsAndCounts.stream().map(pair -> pair.getRight()).collect(Collectors.toList());

        List<double[]> responsibilities = tumorLodsAndCounts.stream()
                .mapToDouble(pair -> MathUtils.arrayMax(pair.getLeft()))
                .mapToObj(lod -> lod < MTFAC.highConfidenceLod ? new double[] {1,0} : new double[] {0,1})
                .collect(Collectors.toList());
        updatePriors(responsibilities);
        fitShape(allCounts, responsibilities);

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
        int g = 90;
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

    public double getSomaticProbability(final double tumorLog10Odds, final double refCount, final double altCount) {
        return 1 - getResponsibilities(tumorLog10Odds, refCount, altCount)[0];
    }

    public boolean passesThreshold(final double tumorLog10Odds, final double refCount, final double altCount) {
        final double nonSomaticProbability =  getResponsibilities(tumorLog10Odds, refCount, altCount)[0];
        return nonSomaticProbability < posteriorThreshold;
    }

    private void updatePriors(final List<double[]> responsibilities) {
        final double somaticCount = responsibilities.stream().mapToDouble(r -> r[1]).sum();
        log10SomaticPrior = FastMath.log10(somaticCount / callableSites);
        log10NothingPrior = FastMath.log10((callableSites - somaticCount)/ callableSites);
    }

    private BetaDistributionShape fitShape(List<double[]> counts, final double[] responsibilities) {
        final int N = counts.size();
        final double weightedAltCount = IntStream.range(0, N).mapToDouble(n -> counts.get(n)[1]*responsibilities[n]).sum();
        final double weightedRefCount = IntStream.range(0, N).mapToDouble(n -> counts.get(n)[0]*responsibilities[n]).sum();
        final double mean = (weightedAltCount + 0.5) / (weightedAltCount + weightedRefCount + 1);
        final double[] totalCounts = IntStream.range(0, N).mapToDouble(n -> MathUtils.sum(counts.get(n))).toArray();

        final double lhs = (1/(mean*(1-mean))) * IntStream.range(0, N).mapToDouble(n ->
                responsibilities[n]*MathUtils.square(counts.get(n)[1] - mean * totalCounts[n])).sum();

        final UnivariateRealFunction lhsMinusRhs = a -> {
                final double b = getBeta(a, mean);
                return lhs - IntStream.range(0,N).mapToDouble(n -> responsibilities[n]*totalCounts[n]*(a+b+totalCounts[n])).sum()/(a+b+1);
        };

        try {
            final double alpha = new BisectionSolver().solve(2000, lhsMinusRhs, 1, 100);
            return new BetaDistributionShape(alpha, getBeta(alpha, mean));
        } catch (MathException ex) {
            throw new GATKException(ex.getMessage());
        }
    }

    private void fitShape(final List<double[]> counts, final List<double[]> responsibilities) {
        Utils.validateArg(counts.size() == responsibilities.size(), "must have one responsibility per count");
        final double[] somaticResponsibilities = responsibilities.stream().mapToDouble(r -> r[1]).toArray();
        somaticDistribution = fitShape(counts, somaticResponsibilities);
    }

    private double getBeta(final double alpha, final double mean) {
        return (1 - mean) * alpha / mean;
    }
}
