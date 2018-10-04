package org.broadinstitute.hellbender.utils;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math3.util.MathArrays;

import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 *
 * @param <CLUSTER>   The type of latent cluster "center" parameters
 * @param <DATUM>     The type of observed data
 */
public abstract class DirichletClusterer<CLUSTER, DATUM> {
    final List<DATUM> data;
    List<CLUSTER> clusters;
    final Dirichlet weightsPrior;
    final RealMatrix responsibilities;

    private static final int NUM_ITERATIONS = 10;

    public DirichletClusterer(final List<DATUM> data, final double concentration) {
        this.data = data;
        clusters = initializeClusters(data);
        weightsPrior = Dirichlet.symmetricDirichlet(clusters.size(), concentration);
        responsibilities = new Array2DRowRealMatrix(data.size(), clusters.size());

        for (int n = 0; n < NUM_ITERATIONS; n++) {
            iterate();
        }

    }

    public Function<DATUM, Double> mixtureModelLogLikelihood() {
        final double[] log10EffectiveWeights = MathUtils.normalizeLog10(weightsPosterior().effectiveLog10MultinomialWeights());
        final double[] logEffectiveWeights = MathUtils.applyToArrayInPlace(log10EffectiveWeights, x -> x * MathUtils.LOG10_OF_E);

        return new Function<DATUM, Double>() {
            final private double[] logWeights = logEffectiveWeights;
            
            public Double apply(final DATUM datum) {
                final double[] logLikelihoods = clusters.stream().mapToDouble(cluster -> logLikelihood(cluster, datum)).toArray();
                return MathUtils.logSumExp(MathArrays.ebeAdd(logWeights, logLikelihoods));
            }
        };
    }

    private void iterate() {
        relearnResponsibilities();
        relearnClusters();
    }

    private void relearnResponsibilities() {
        final double[] logEffectiveWeights = weightsPosterior().effectiveLogMultinomialWeights();
        for (int datumIndex = 0; datumIndex < data.size(); datumIndex++) {
            final DATUM datum = data.get(datumIndex);
            final double[] logLikelihoods = clusters.stream().mapToDouble(cluster -> logLikelihood(cluster, datum)).toArray();
            final double[] log10Posteriors = MathUtils.applyToArrayInPlace(MathArrays.ebeAdd(logEffectiveWeights, logLikelihoods), x -> x * MathUtils.LOG10_OF_E);
            final double[] datumRepsonsibilities = MathUtils.normalizeFromLog10ToLinearSpace(log10Posteriors);
            responsibilities.setRow(datumIndex, datumRepsonsibilities);
        }
    }

    private void relearnClusters() {
        clusters = IntStream.range(0, clusters.size()).mapToObj(k -> relearnCluster(clusters.get(k), data, responsibilities.getColumn(k)))
                .collect(Collectors.toList());

    }

    public Dirichlet weightsPosterior() {
        final double[] pseudocounts =  new IndexRange(0, clusters.size()).mapToDouble(k -> MathUtils.sum(responsibilities.getColumn(k)));
        return weightsPrior.addCounts(pseudocounts);
    }

    public abstract List<CLUSTER> initializeClusters(final List<DATUM> data);

    public abstract double logLikelihood(final CLUSTER cluster, final DATUM datum);

    public abstract CLUSTER relearnCluster(final CLUSTER current, final List<DATUM> data, final double[] responsibilities);
}
