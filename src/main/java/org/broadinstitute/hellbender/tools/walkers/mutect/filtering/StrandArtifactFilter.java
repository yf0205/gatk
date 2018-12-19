package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import com.google.common.primitives.Doubles;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandArtifact;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Arrays;
import java.util.List;

// TODO: this should not be a hard filter because it emits a posterior!!!
public class StrandArtifactFilter extends HardFilter {
    @Override
    public boolean isArtifact(final VariantContext vc, final Mutect2FilteringInfo filteringInfo) {
        final List<Double> posteriorProbabilities = vc.getAttributeAsDoubleList(GATKVCFConstants.STRAND_ARTIFACT_POSTERIOR_KEY, 0.0);
        final List<Double> mapAlleleFractionEstimates = vc.getAttributeAsDoubleList(GATKVCFConstants.STRAND_ARTIFACT_AF_KEY, 0.0);

        final int maxZIndex = MathUtils.maxElementIndex(Doubles.toArray(posteriorProbabilities));

        if (maxZIndex == StrandArtifact.ArtifactState.NO_ARTIFACT.ordinal()){
            return false;
        }

        return posteriorProbabilities.get(maxZIndex) > filteringInfo.getMTFAC().strandArtifactPosteriorProbThreshold &&
                mapAlleleFractionEstimates.get(maxZIndex) < filteringInfo.getMTFAC().strandArtifactAlleleFractionThreshold;
    }

    public String filterName() {
        return GATKVCFConstants.STRAND_ARTIFACT_FILTER_NAME;
    }

    protected List<String> requiredAnnotations() {
        return Arrays.asList(GATKVCFConstants.STRAND_ARTIFACT_POSTERIOR_KEY, GATKVCFConstants.STRAND_ARTIFACT_AF_KEY);
    }
}
