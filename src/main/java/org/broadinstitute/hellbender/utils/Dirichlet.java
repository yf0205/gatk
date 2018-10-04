package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.Arrays;
import java.util.Collections;
import java.util.stream.IntStream;

/**
 * The Dirichlet distribution is a distribution on multinomial distributions: if pi is a vector of positive multinomial weights
 * such that sum_i pi[i] = 1, the Dirichlet pdf is P(pi) = [prod_i Gamma(alpha[i]) / Gamma(sum_i alpha[i])] * prod_i pi[i]^(alpha[i] - 1)
 *
 * The vector alpha comprises the sufficient statistics for the Dirichlet distribution.
 *
 * Since the Dirichlet is the conjugate prior to the multinomial, if one has a Dirichlet prior with concentration alpha
 * and observes each category i n_i times (assuming categories are drawn from a multinomial distribution pi)
 * the posterior is alpha_i -> alpha_i + n_i
 *
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public class Dirichlet {
    final protected double[] alpha;

    protected Dirichlet(final double... alpha) {
        Utils.nonNull(alpha);
        Utils.validateArg(alpha.length >= 1, "Dirichlet parameters must have at least one element");
        Utils.validateArg(MathUtils.allMatch(alpha, x -> x >= 0), "Dirichlet parameters may not be negative");
        Utils.validateArg(MathUtils.allMatch(alpha, Double::isFinite), "Dirichlet parameters must be finite");
        this.alpha = alpha.clone();
    }

    public static Dirichlet of(final double... alpha) {
        return new Dirichlet(alpha);
    }

    public static Dirichlet flat(final int dimension) {
        ParamUtils.isPositive(dimension, "Dimension must be positive");
        return of(new IndexRange(0, dimension).mapToDouble(n -> 1));
    }

    /**
     * Create a symmetric distribution Dir(a/K, a/K, a/K . . .) where K is the number of states and
     * a is the concentration.
     */
    public static Dirichlet symmetricDirichlet(final int numStates, final double concentration) {
        Utils.validateArg(numStates > 0, "Must have at leat one state");
        Utils.validateArg(concentration > 0, "concentration must be positive");
        return of(Collections.nCopies(numStates, concentration/numStates).stream().mapToDouble(x->x).toArray());
    }

    // in variational Bayes one often needs the effective point estimate of a multinomial distribution with a
    // Dirichlet prior.  This value is not the mode or mean of the Dirichlet but rather the exp of the expected log weights.
    // note that these effective weights do not add up to 1.  This is fine because in any probabilistic model scaling all weights
    // amounts to an arbitrary normalization constant, but it's important to keep in mind because some classes may expect
    // normalized weights.  In that case the calling code must normalize the weights.
    public double[] effectiveMultinomialWeights() {
        final double digammaOfSum = Gamma.digamma(MathUtils.sum(alpha));
        return MathUtils.applyToArray(alpha, a -> Math.exp(Gamma.digamma(a) - digammaOfSum));
    }

    public double[] effectiveLog10MultinomialWeights() {
        final double digammaOfSum = Gamma.digamma(MathUtils.sum(alpha));
        return MathUtils.applyToArray(alpha, a -> (Gamma.digamma(a) - digammaOfSum) * MathUtils.LOG10_OF_E);
    }

    public double[] effectiveLogMultinomialWeights() {
        final double digammaOfSum = Gamma.digamma(MathUtils.sum(alpha));
        return MathUtils.applyToArray(alpha, a -> Gamma.digamma(a) - digammaOfSum);
    }

    // the array of E[ln x_i] for components i
    public double[] meanLogs() {
        final double digammaOfSum = Gamma.digamma(MathUtils.sum(alpha));
        return MathUtils.applyToArray(alpha, a -> Gamma.digamma(a) - digammaOfSum);
    }

    public double[] meanWeights() {
        final double sum = MathUtils.sum(alpha);
        return MathUtils.applyToArray(alpha, x -> x / sum);
    }

    public double[] log10MeanWeights() {
        final double sum = MathUtils.sum(alpha);
        return MathUtils.applyToArray(alpha, x -> Math.log10(x / sum));
    }

    public double logNormalization() {
        final double logNumerator = Gamma.logGamma(MathUtils.sum(alpha));
        final double logDenominator = MathUtils.sum(MathUtils.applyToArray(alpha, Gamma::logGamma));
        return logNumerator - logDenominator;
    }

    public double log10Normalization() {
        return MathUtils.logToLog10(logNormalization());
    }

    public int dimension() { return alpha.length; }

    public double distance1(final Dirichlet other) {
        Utils.validateArg(this.dimension() == other.dimension(), "Dirichlets must have same dimension.");
        return MathArrays.distance1(this.alpha, other.alpha);
    }

    public Dirichlet addCounts(final double[] counts) {
        Utils.validateArg(counts.length == dimension(), "Counts must have same dimension as prior.");
        return of(MathArrays.ebeAdd(alpha, counts));
    }
}
