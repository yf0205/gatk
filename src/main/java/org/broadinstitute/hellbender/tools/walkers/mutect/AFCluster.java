package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.BetaBinomialDistribution;
import org.broadinstitute.hellbender.utils.BetaDistributionShape;
import org.broadinstitute.hellbender.utils.Dirichlet;

public class AFCluster {
    final private boolean isBackground;

    final private BetaDistributionShape shape;

    public AFCluster(final boolean isBackground, final BetaDistributionShape shape) {
        this.isBackground = isBackground;
        this.shape = shape;
    }

    public boolean isBackground() { return isBackground; }

    public double logLikelihood(final Count count) {
        return new BetaBinomialDistribution(null, shape.getAlpha(), shape.getBeta(), count.getTotal()).logProbability(count.getAltCount());
    }

    public Dirichlet asDirichlet() { return Dirichlet.of(shape.getAlpha(), shape.getBeta()); }
}
