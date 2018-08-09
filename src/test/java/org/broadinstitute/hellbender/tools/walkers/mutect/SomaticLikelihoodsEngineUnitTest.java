package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.utils.Dirichlet;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.stream.IntStream;

/**
 * Created by David Benjamin on 3/9/17.
 */
public class SomaticLikelihoodsEngineUnitTest extends GATKBaseTest {

    @DataProvider(name="afPosterior")
    public Object[][] afPosterior() {
        // Dirichlet prior, log 10 likelihoods matrix, expected counts, whether to test log 10 evidence correction
        return new Object[][]{
                //likelihoods completely favor allele 0 over allele 1 for every read, so we should get no counts for allele 1
                { Dirichlet.of(1, 1),
                        new Array2DRowRealMatrix(new double[][] {{0, 0, 0, 0}, {-10, -10, -10, -10}}),
                        new double[] {4, 0},
                        true},
                //extremely strong prior, extremely weak likelihoods
                { Dirichlet.of(1e8, 1),
                        new Array2DRowRealMatrix(new double[][] {{0, 0, 0, 0}, {0, 0, 0, 0}}),
                        new double[] {4, 0},
                        false},
                // extremely weak prior, extremely strong likelihoods
                { Dirichlet.of(1e-6, 1e-6),
                        new Array2DRowRealMatrix(new double[][] {{0, 0, 0, -10}, {-10, -10, -10, 0}}),
                        new double[] {3, 1},
                        true},
                // non-obvious expected counts, but we can still test convergence
                { Dirichlet.of(0.2, 1.7),
                        new Array2DRowRealMatrix(new double[][] {{0.1, 5.2, 0.5, 0.2}, {2.6, 0.6, 0.5, 0.4}}),
                        null,
                        true},
        };
    }

    @Test(dataProvider = "afPosterior")
    public void testAlleleFractionsPosterior(final Dirichlet prior, final RealMatrix matrix, final double[] expectedCounts, final boolean testEvidenceCorrection) {
        final Dirichlet posterior = SomaticLikelihoodsEngine.alleleFractionsPosterior(matrix, prior);

        if (expectedCounts != null) {
            final Dirichlet expectedPosterior = prior.addCounts(expectedCounts);
            Assert.assertEquals(posterior.distance1(expectedPosterior), 0, 1.0e-6);
        }

        // test convergence i.e. posterior = prior + effective counts
        final double[] effectiveCounts = SomaticLikelihoodsEngine.getEffectiveCounts(matrix, posterior);
        Assert.assertEquals(prior.addCounts(effectiveCounts).distance1(posterior), 0, 1.0e-3);
    }

    @Test(dataProvider = "afPosterior")
    public void testLog10EvidenceCorrection(final Dirichlet prior, final RealMatrix log10Likelihoods, final double[] expectedCounts, final boolean testEvidenceCorrection) {
        if (!testEvidenceCorrection) {
            return;
        }
        final double evidence = SomaticLikelihoodsEngine.log10Evidence(log10Likelihoods, prior);
        final Dirichlet posterior = SomaticLikelihoodsEngine.alleleFractionsPosterior(log10Likelihoods, prior);
        final double[] counts = SomaticLikelihoodsEngine.getEffectiveCounts(log10Likelihoods, posterior);

        final Dirichlet otherPrior = Dirichlet.of(IntStream.range(0, prior.dimension()).mapToDouble(n -> n + 3).toArray());
        final double otherEvidence = SomaticLikelihoodsEngine.log10Evidence(log10Likelihoods, otherPrior);

        final double correction = SomaticLikelihoodsEngine.log10OddsCorrection(otherPrior, prior, counts);

        Assert.assertEquals(evidence + correction, otherEvidence, 0.1);
    }

    @Test
    public void testEvidence() {
        // one exact limit for the evidence is when the likelihoods of each read are so peaked (i.e. the most likely allele
        // of each read is much likelier than all other alleles) that the sum over latent read-to-allele assignments
        // (that is, over the indicator z in the notes) is dominated by the max-likelihood allele configuration
        // and thus the evidence reduces to exactly integrating out the Dirichlet allele fractions

        final Dirichlet prior = Dirichlet.of(1,2);
        final RealMatrix log10Likelihoods = new Array2DRowRealMatrix(new double[][] {{0.1, 4.0, 3.0, -10}, {-12, -9, -5.0, 0.5}});
        final double calculatedLog10Evidence = SomaticLikelihoodsEngine.log10Evidence(log10Likelihoods, prior);
        final double[] maxLikelihoodCounts = new double[] {3, 1};
        final double expectedLog10Evidence = prior.log10Normalization() - prior.addCounts(maxLikelihoodCounts).log10Normalization()
                + new IndexRange(0,log10Likelihoods.getColumnDimension()).sum(read -> log10Likelihoods.getColumnVector(read).getMaxValue());
        Assert.assertEquals(calculatedLog10Evidence, expectedLog10Evidence, 1e-5);

        // when there's just one read we can calculate the likelihood exactly

        final Dirichlet prior2 = Dirichlet.of(1,2);
        final RealMatrix log10Likelihoods2 = new Array2DRowRealMatrix(new double[][] {{0.1}, {0.5}});
        final double calculatedLog10Evidence2 = SomaticLikelihoodsEngine.log10Evidence(log10Likelihoods2, prior2);
        final double[] delta0 = new double[] {1, 0};
        final double[] delta1 = new double[] {0, 1};
        final double expectedLog10Evidence2 = MathUtils.log10SumLog10(log10Likelihoods2.getEntry(0,0) +
                prior2.log10Normalization() - prior2.addCounts(delta0).log10Normalization(),
                log10Likelihoods2.getEntry(1,0) +
                prior2.log10Normalization() - prior2.addCounts(delta1).log10Normalization());
        Assert.assertEquals(calculatedLog10Evidence2, expectedLog10Evidence2, 0.05);


    }

}