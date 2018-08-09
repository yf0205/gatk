package org.broadinstitute.hellbender.utils;

import org.broadinstitute.hellbender.utils.Dirichlet;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

public class BetaDistributionShape extends Dirichlet {

    public BetaDistributionShape(final double alpha, final double beta){
        super(alpha, beta);
    }

    public double getAlpha() { return alpha[0]; }

    public double getBeta() {
        return alpha[1];
    }
}
