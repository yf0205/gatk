package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Arrays;
import java.util.List;

public class ShortTandemRepeatFilter extends HardFilter {
    @Override
    public boolean isArtifact(final VariantContext vc, final Mutect2FilteringInfo filteringInfo) {

        final int[] rpa = vc.getAttributeAsList(GATKVCFConstants.REPEATS_PER_ALLELE_KEY).stream()
                .mapToInt(o -> Integer.parseInt(String.valueOf(o))).toArray();
        if (rpa.length < 2) {
            return false;
        }
        final String ru = vc.getAttributeAsString(GATKVCFConstants.REPEAT_UNIT_KEY, "");

        final int referenceSTRBaseCount = ru.length() * rpa[0];
        final int numPCRSlips = rpa[0] - rpa[1];
        if (referenceSTRBaseCount >= filteringInfo.getMTFAC().minPcrSlippageBases && Math.abs(numPCRSlips) == 1) {
            // calculate the p-value that out of n reads we would have at least k slippage reads
            // if this p-value is small we keep the variant (reject the PCR slippage hypothesis)
            final int[] ADs = sumADsOverSamples(vc, filteringInfo.getNormalSamples(), true, false);

            final int depth = ADs == null ? 0 : (int) MathUtils.sum(ADs);
            final double oneSidedPValueOfSlippage = (ADs == null || ADs.length < 2) ? 1.0 :
                    new BinomialDistribution(null, depth, filteringInfo.getMTFAC().pcrSlippageRate).cumulativeProbability(ADs[1] - 1, depth);
            return oneSidedPValueOfSlippage > filteringInfo.getMTFAC().pcrSlippagePValueThreshold;
        } else {
            return false;
        }
    }

    public String filterName() {
        return GATKVCFConstants.STR_CONTRACTION_FILTER_NAME;
    }

    protected List<String> requiredAnnotations() {
        return Arrays.asList(GATKVCFConstants.REPEATS_PER_ALLELE_KEY, GATKVCFConstants.REPEAT_UNIT_KEY);
    }
}
