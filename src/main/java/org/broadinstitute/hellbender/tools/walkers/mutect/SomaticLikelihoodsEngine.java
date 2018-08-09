package org.broadinstitute.hellbender.tools.walkers.mutect;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.utils.*;

import java.util.Arrays;

/**
 * Created by David Benjamin on 3/9/17.
 */
public class SomaticLikelihoodsEngine {

    public static final double CONVERGENCE_THRESHOLD = 0.001;

    /**
     * Given a likelihoods matrix, calculate the parameters of the Dirichlet posterior distribution on their allele
     * fractions, which define a discrete distribution.
     * @param log10Likelihoods matrix of alleles x reads
     * @param prior
     */
    @VisibleForTesting
    static Dirichlet alleleFractionsPosterior(final RealMatrix log10Likelihoods, final Dirichlet prior) {
        final int numberOfAlleles = log10Likelihoods.getRowDimension();
        Utils.validateArg(numberOfAlleles == prior.dimension(), "Must have one pseudocount per allele.");

        Dirichlet posterior = Dirichlet.flat(numberOfAlleles);
        boolean converged = false;

        while(!converged) {
            // alleleCounts = \sum_r \bar{z}_r, where \bar{z}_r is an a-dimensional vector of the expectation of z_r with respect to q(f)
            final double[] alleleCounts = getEffectiveCounts(log10Likelihoods, posterior);
            final Dirichlet newPosterior = prior.addCounts(alleleCounts);
            converged = newPosterior.distance1(posterior) < CONVERGENCE_THRESHOLD;
            posterior = newPosterior;
        }

        return posterior;
    }

    /**
     * Given data log likelihoods and a Dirichlet prior for a categorical distribution, obtain the array of total
     * responsibilities for each category
     * @param log10Likelihoods
     * @param prior Dirichlet prior on allele fractions
     * @return
     */
    @VisibleForTesting
    protected static double[] getEffectiveCounts(RealMatrix log10Likelihoods, final Dirichlet prior) {
        final double[] effectiveLog10Weights = prior.effectiveLog10MultinomialWeights();
        return MathUtils.sumArrayFunction(0, log10Likelihoods.getColumnDimension(),
                read -> MathUtils.posteriors(effectiveLog10Weights, log10Likelihoods.getColumn(read)));
    }


    /**
     * @param log10Likelihoods matrix of alleles x reads
     * @param prior Dirichlet prior on allele fractions
     */
    @VisibleForTesting
    static double log10Evidence(final RealMatrix log10Likelihoods, final Dirichlet prior) {
        Utils.validateArg(log10Likelihoods.getRowDimension() == prior.dimension(), "Must have one pseudocount per allele.");
        final Dirichlet posterior = alleleFractionsPosterior(log10Likelihoods, prior);

        final double[] log10AlleleFractions = posterior.effectiveLog10MultinomialWeights();
        final double likelihoodsAndEntropyContribution = new IndexRange(0, log10Likelihoods.getColumnDimension()).sum(r -> {
            final double[] log10LikelihoodsForRead = log10Likelihoods.getColumn(r);
            final double[] responsibilities = MathUtils.posteriors(log10AlleleFractions, log10LikelihoodsForRead);
            final double likelihoodsTerm = MathUtils.sum(MathArrays.ebeMultiply(log10LikelihoodsForRead, responsibilities));
            final double entropyTerm = Arrays.stream(responsibilities).map(MathUtils::xLog10x).sum();
            return likelihoodsTerm - entropyTerm;

        });

        return prior.log10Normalization() + -posterior.log10Normalization() + likelihoodsAndEntropyContribution;
    }

    // same as above using the default flat prior
    public static double log10Evidence(final RealMatrix log10Likelihoods) {
        return log10Evidence(log10Likelihoods, Dirichlet.flat(log10Likelihoods.getRowDimension()));
    }

    // correct the log 10 odds originally output by Mutect2, which assumed a flat prior, to approximate
    // the log 10 odds one would obtain from using a different, non-flat prior
    public static double log10OddsCorrection(final Dirichlet newPrior, final Dirichlet oldPrior, final double[] counts) {
        Utils.validateArg(newPrior.dimension() == oldPrior.dimension(), "Priors must have same dimensions");
        Utils.validateArg(newPrior.dimension() == counts.length, "Priors must have same dimensions as counts.");

        final double oldValue = oldPrior.log10Normalization() - oldPrior.addCounts(counts).log10Normalization();
        final double newValue = newPrior.log10Normalization() - newPrior.addCounts(counts).log10Normalization();
        return newValue - oldValue;
    }
}
