package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;

public class ChimericOriginalAlignmentFilter extends HardFilter {
    @Override
    public boolean isArtifact(final VariantContext vc, final Mutect2FilteringInfo filteringInfo) {
        if(!vc.isBiallelic()) {
            return false;
        }

        final int altCount = vc.getGenotypes().stream().mapToInt(g -> g.getAD()[1]).sum();
        final int nonMtOa = vc.getAttributeAsInt(GATKVCFConstants.ORIGINAL_CONTIG_MISMATCH_KEY, 0);
        return (double) nonMtOa / altCount > filteringInfo.getMTFAC().nonMtAltByAlt;
    }

    public String filterName() {
        return GATKVCFConstants.CHIMERIC_ORIGINAL_ALIGNMENT_FILTER_NAME;
    }

    protected List<String> requiredAnnotations() { return Collections.singletonList(GATKVCFConstants.ORIGINAL_CONTIG_MISMATCH_KEY); }
}
