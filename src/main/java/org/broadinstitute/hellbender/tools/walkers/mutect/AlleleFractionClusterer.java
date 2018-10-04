package org.broadinstitute.hellbender.tools.walkers.mutect;

import com.google.common.collect.Lists;
import org.apache.commons.math.MathException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.solvers.BisectionSolver;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.*;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class AlleleFractionClusterer extends DirichletClusterer<AFCluster, Count> {
    private static final double MIN_SHAPE_SIZE_FOR_BINOMIAL = 1000;

    private static final int NUM_CLUSTERS = 20;

    public AlleleFractionClusterer(final List<Count> data, final double concentration) {
        super(data, concentration);
    }

    public double log10OddsCorrection(final double altCount, final double refCount) {
        final double[] log10ClusterCorrections = clusters.stream()
                .mapToDouble(cluster -> SomaticLikelihoodsEngine.log10OddsCorrection(
                cluster.asDirichlet(), Dirichlet.flat(2), new double[]{altCount, refCount}))
                .toArray();

        return MathUtils.log10SumLog10(MathArrays.ebeAdd(log10ClusterCorrections, log10EffectiveWeights));
    }

    public double logLikelihood(final AFCluster cluster, final Count count) {
        return cluster.logLikelihood(count);
    }

    public List<AFCluster> initializeClusters(final List<Count> data) {
        //sort by AF
        final List<Double> sortedAFs = data.stream().map(Count::getAF).sorted().collect(Collectors.toList());

        final int chunkSize = Math.max(sortedAFs.size() / NUM_CLUSTERS, 1);
        final List<Double> chunkAverageAFs = Lists.partition(sortedAFs, chunkSize).stream()
                .map(chunk -> chunk.stream().mapToDouble(af -> af).average().getAsDouble())
                .distinct()
                .collect(Collectors.toList());

        // add an initially flat background Beta cluster last
        final List<AFCluster> clusters = chunkAverageAFs.stream().map(af -> makeBinomialCluster(af)).collect(Collectors.toList());
        clusters.add(new AFCluster(true, new BetaDistributionShape(1,1)));
        return clusters;
    }

    public AFCluster relearnCluster(final AFCluster current, final List<Count> data, final double[] responsibilities) {
        if (current.isBackground()) {
            return new AFCluster(true, fitShape(data, responsibilities));
        } else {    // weighted average of counts
            final double weightedAltCount = new IndexRange(0, data.size()).sum(n -> responsibilities[n] * data.get(n).getAltCount());
            final double weightedRefCount = new IndexRange(0, data.size()).sum(n -> responsibilities[n] * data.get(n).getRefCount());
            final double alleleFraction = weightedAltCount / (weightedAltCount + weightedRefCount);
            return makeBinomialCluster(alleleFraction);
        }
    }

    private AFCluster makeBinomialCluster(final double alleleFraction) {
        final double alpha = MIN_SHAPE_SIZE_FOR_BINOMIAL;
        final double beta = alpha * (1 - alleleFraction) / alleleFraction;
        return new AFCluster(false, new BetaDistributionShape(alpha, beta));
    }

    private BetaDistributionShape fitShape(List<Count> counts, final double[] responsibilities) {
        final int N = counts.size();
        final double weightedAltCount = new IndexRange(0, N).sum(n -> counts.get(n).getAltCount()*responsibilities[n]);
        final double weightedRefCount = new IndexRange(0, N).sum(n -> counts.get(n).getRefCount()*responsibilities[n]);
        final double mean = (weightedAltCount + 0.5) / (weightedAltCount + weightedRefCount + 1);
        final double[] totalCounts = counts.stream().mapToDouble(Count::getTotal).toArray();

        final double lhs = (1/(mean*(1-mean))) * new IndexRange(0, N).sum(n ->
                responsibilities[n]*MathUtils.square(counts.get(n).getRefCount() - mean * totalCounts[n]));

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

    private double getBeta(final double alpha, final double mean) {
        return (1 - mean) * alpha / mean;
    }

}
