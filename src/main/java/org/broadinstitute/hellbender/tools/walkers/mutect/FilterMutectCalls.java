package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import it.unimi.dsi.fastutil.Hash;
import org.apache.commons.lang.mutable.MutableInt;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.exome.FilterByOrientationBias;
import org.broadinstitute.hellbender.tools.walkers.contamination.CalculateContamination;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * <p>Filter variants in a Mutect2 VCF callset.</p>
 *
 * <p>
 *     FilterMutectCalls encapsulates GATK3 MuTect2's filtering functionality and adds additional filters.
 *     Thresholds for filters are contained in {@link M2FiltersArgumentCollection} and described in
 *     <a href='https://github.com/broadinstitute/gatk/tree/master/docs/mutect/mutect.pdf' target='_blank'>https://github.com/broadinstitute/gatk/tree/master/docs/mutect/mutect.pdf</a>.
 *     To filter based on sequence context artifacts, see {@link FilterByOrientationBias}.
 * </p>
 * <p>
 *     Filtering thresholds for both normal-artifact-lod (default threshold 0.0) and tumor-lod (default threshold 5.3) can be set in this tool.
 *     If the normal artifact log odds is larger than the threshold, then FilterMutectCalls applies the artifact-in-normal filter.
 *     For matched normal analyses with tumor contamination in the normal, consider increasing the normal-artifact-lod threshold.
 *     If the tumor log odds is smaller than the threshold, then FilterMutectCalls filters the variant.
 * </p>
 * <p>
 *     If given a --contamination-table file, e.g. results from
 *     {@link CalculateContamination}, the tool will additionally
 *     filter on contamination fractions. Alternatively, provide a numerical fraction to filter with the --contamination argument.
 * </p>
 * <p>
 *     This tool is featured in the Somatic Short Mutation calling Best Practice Workflow.
 *     See <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11136">Tutorial#11136</a> for a
 *     step-by-step description of the workflow and <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11127">Article#11127</a>
 *     for an overview of what traditional somatic calling entails. For the latest pipeline scripts, see the
 *     <a href="https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl">Mutect2 WDL scripts directory</a>.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk FilterMutectCalls \
 *   -V somatic.vcf.gz \
 *   --contamination-table contamination.table \
 *   -O filtered.vcf.gz
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Filter somatic SNVs and indels called by Mutect2",
        oneLineSummary = "Filter somatic SNVs and indels called by Mutect2",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
public final class FilterMutectCalls extends TwoPassVariantWalker {

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output filtered VCF file", optional=false)
    private final String outputVcf = null;

    @ArgumentCollection
    protected M2FiltersArgumentCollection MTFAC = new M2FiltersArgumentCollection();

    private VariantContextWriter vcfWriter;

    private Map<String, Pair<Mutect2FilteringEngine, FilteringFirstPass>> filteringBySample;

    @Override
    public void onTraversalStart() {
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> headerLines = inputHeader.getMetaDataInSortedOrder().stream()
                .filter(line -> !line.getKey().equals(Mutect2FilteringEngine.FILTERING_STATUS_VCF_KEY)) //remove header line from Mutect2 stating that calls are unfiltered.
                .collect(Collectors.toSet());
        headerLines.add(new VCFHeaderLine(Mutect2FilteringEngine.FILTERING_STATUS_VCF_KEY, "These calls have been filtered by " + FilterMutectCalls.class.getSimpleName() + " to label false positives with a list of failed filters and true positives with PASS."));

        GATKVCFConstants.MUTECT_FILTER_NAMES.stream().map(GATKVCFHeaderLines::getFilterLine).forEach(headerLines::add);

        headerLines.addAll(getDefaultToolVCFHeaderLines());

        final VCFHeader vcfHeader = new VCFHeader(headerLines, inputHeader.getGenotypeSamples());
        vcfWriter = createVCFWriter(new File(outputVcf));
        vcfWriter.writeHeader(vcfHeader);

        final Optional<String> normalSample = getNormalSampleName();

        final List<String> tumorSamples = getHeaderForVariants().getGenotypeSamples().stream()
                .filter(sample -> !(normalSample.isPresent() && sample.equals(normalSample.get())))
                .collect(Collectors.toList());

        filteringBySample = tumorSamples.stream()
                .collect(Collectors.toMap(sample -> sample, sample -> ImmutablePair.of(new Mutect2FilteringEngine(MTFAC, sample, normalSample), new FilteringFirstPass(sample))));
    }

    private Optional<String> getNormalSampleName() {
        final VCFHeaderLine normalSampleHeaderLine = getHeaderForVariants().getMetaDataLine(Mutect2Engine.NORMAL_SAMPLE_KEY_IN_VCF_HEADER);
        return normalSampleHeaderLine == null ? Optional.empty() : Optional.of(normalSampleHeaderLine.getValue());
    }

    @Override
    public Object onTraversalSuccess() {
        return "SUCCESS";
    }

    @Override
    public void firstPassApply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext refContext, final FeatureContext fc) {
        filteringBySample.values().forEach(pair -> pair.getRight().add(pair.getLeft().calculateFilters(MTFAC, vc, Optional.empty()), vc));
    }

    @Override
    protected void afterFirstPass() {
        final List<FilteringFirstPass> firstPasses = filteringBySample.values().stream().map(pair -> pair.getRight()).collect(Collectors.toList());
        firstPasses.forEach(firstPass -> firstPass.learnModelForSecondPass(MTFAC.maxFalsePositiveRate));
        FilteringFirstPass.writeM2FilterSummary(firstPasses, MTFAC.mutect2FilteringStatsTable);
    }

    @Override
    public void secondPassApply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext refContext, final FeatureContext fc) {
        final MutableInt totalAltDepth = new MutableInt(0);
        final Map<String, MutableInt> filteredCountsByFilter = new HashMap<>();
        final MutableInt filteredCounts = new MutableInt(0);

        for (final Map.Entry<String, Pair<Mutect2FilteringEngine, FilteringFirstPass>> entry : filteringBySample.entrySet()) {
            final String sample = entry.getKey();
            final Mutect2FilteringEngine engine = entry.getValue().getLeft();
            final FilteringFirstPass firstPass = entry.getValue().getRight();

            final int[] alleleDepths = vc.getGenotype(sample).getAD();
            final int altDepth = (int) MathUtils.sum(alleleDepths) - alleleDepths[0];

            if (altDepth == 0) {
                continue;
            }

            final FilterResult result = engine.calculateFilters(MTFAC, vc, Optional.of(firstPass));

            totalAltDepth.add(altDepth);
            if (!result.getFilters().isEmpty()) {
                filteredCounts.add(altDepth);
            }
            result.getFilters().forEach(filter -> {
                filteredCountsByFilter.putIfAbsent(filter, new MutableInt(0));
               filteredCountsByFilter.get(filter).add(altDepth);
            });
        }

        final VariantContextBuilder vcb = new VariantContextBuilder(vc).filters(new HashSet<>());;
        if (filteredCounts.intValue() > MTFAC.filteredSampleFraction * totalAltDepth.intValue()) {
            vcb.filters(filteredCountsByFilter.keySet());
        }

        vcfWriter.add(vcb.make());
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
