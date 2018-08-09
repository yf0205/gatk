package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.broadinstitute.hellbender.tools.walkers.mutect.SomaticLikelihoodsEngine;
import org.testng.Assert;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class DirichletUnitTest {

    @Test
    public void testDirichletNormalization() {
        // a Dirichlet with two parameters is a Beta
        final double a = 11;
        final double b = 36;
        final Dirichlet params1 = Dirichlet.of(a, b);
        final double normalization = Math.pow(10, params1.log10Normalization());
        final BetaDistribution bd = new BetaDistribution(a,b);
        for (final double x : new double[] {0.1, 0.3, 0.6}) {
            Assert.assertEquals(bd.density(x), normalization * Math.pow(x, a - 1) * Math.pow(1-x, b-1), 1e-6);
        }

        // a Dirichlet with parameters equal to 1 is flat, so the normalization is 1/the hypervolume of the simplex
        // which is d! where d is the dimension

        final Dirichlet params2 = Dirichlet.of(1, 1, 1, 1);
        final double normalization2 = Math.pow(10, params2.log10Normalization());
        Assert.assertEquals(normalization2, 6, 1e-6);
    }

    @Test
    public void testNormalization() {
        final double alpha = 5.7;
        final double beta = 9.9;
        final double x = 0.3;

        final double normalization = Math.pow(10,Dirichlet.of(alpha, beta).log10Normalization());
        final double rest = Math.pow(x, alpha - 1) * Math.pow(1 - x, beta - 1);
        Assert.assertEquals(normalization * rest, new BetaDistribution(null, alpha, beta).density(x), 1e-6);
    }
}