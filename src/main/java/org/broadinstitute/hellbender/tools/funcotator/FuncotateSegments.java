package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;
import org.broadinstitute.hellbender.tools.funcotator.metadata.VcfFuncotationMetadata;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.File;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

// TODO: Fill in the docs and do not forget Oncotator citation.
@CommandLineProgramProperties(
        summary = " (similar functionality to Oncotator).",
        oneLineSummary = "Functional Annotator",
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class FuncotateSegments extends FeatureWalker<AnnotatedInterval> {
    private static final Logger logger = LogManager.getLogger(FuncotateSegments.class);
    @Argument(
            doc = "Input segment file (tab-separated values).",
            fullName = CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME
    )
    private File segmentFile;

    @ArgumentCollection
    private final FuncotatorSegmentArgumentCollection funcotatorArgs = new FuncotatorSegmentArgumentCollection();


    private FuncotatorEngine funcotatorEngine;

    @Override
    protected boolean isAcceptableFeatureType(final Class<? extends Feature> featureType) {
        return featureType.equals(AnnotatedInterval.class);
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    private OutputRenderer outputRenderer;

    @Override
    public void onTraversalStart() {

        Utils.validateArg(funcotatorArgs.transcriptSelectionMode != TranscriptSelectionMode.ALL, "Cannot funcotate segments with the ALL transcript selection mode.  Please select another mode.");

        // Get our overrides for annotations:
        final LinkedHashMap<String, String> annotationDefaultsMap = FuncotatorEngine.splitAnnotationArgsIntoMap(funcotatorArgs.annotationDefaults);
        final LinkedHashMap<String, String> annotationOverridesMap = FuncotatorEngine.splitAnnotationArgsIntoMap(funcotatorArgs.annotationOverrides);

        // Next set up our transcript list:
        final Set<String> finalUserTranscriptIdSet = FuncotatorEngine.processTranscriptList(funcotatorArgs.userTranscriptIdSet);

        // Initialize the datasources (and make sure to filter to handle only segment-enabled funcotation factories.
        // Initialize all of our data sources:
        // Sort data sources to make them process in the same order each time:
        funcotatorArgs.dataSourceDirectories.sort(Comparator.naturalOrder());
        final Map<Path, Properties> configData = DataSourceUtils.getAndValidateDataSourcesFromPaths(funcotatorArgs.referenceVersion, funcotatorArgs.dataSourceDirectories);

        // Create the data sources from the input:
        // This will also create and register the FeatureInputs (created by the Data Sources)
        // with the GATK Engine, so we do not have to plumb them in after the fact.
        //  Only take datasources that support the annotation of segments.
        final List<DataSourceFuncotationFactory> dataSourceFuncotationFactories = DataSourceUtils.createDataSourceFuncotationFactoriesForDataSources(
                configData,
                annotationOverridesMap,
                funcotatorArgs.transcriptSelectionMode,
                finalUserTranscriptIdSet,
                this,
                funcotatorArgs.lookaheadFeatureCachingInBp, new FlankSettings(0,0))
                .stream()
                .filter(ff -> ff.isSupportingSegmentFuncotation())
                .collect(Collectors.toList());

        // Initialize a funcotator engine to handle segments.
        funcotatorEngine = new FuncotatorEngine(funcotatorArgs,
                getBestAvailableSequenceDictionary(),
                VcfFuncotationMetadata.create(
                        Arrays.asList(
                                new VCFInfoHeaderLine(VCFConstants.END_KEY, 1,
                                        VCFHeaderLineType.Integer, "End coordinate of the variant")
                        )
                ),
                dataSourceFuncotationFactories
        );

        // Create a composite output renderer.
        // TODO: Make composite.  For now just using SimpleTsv
        outputRenderer = funcotatorEngine.createOutputRenderer(
                annotationDefaultsMap,
                annotationOverridesMap,
                new VCFHeader(),
                getDefaultToolVCFHeaderLines(),
                this
        );

    }

    @Override
    public void apply(final AnnotatedInterval feature, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        // Convert feature into a VariantContext while honoring the funcotation engine's necessary conversions.
        final VariantContext segmentVariantContext = funcotatorEngine.getDefaultVariantTransformer().apply(
                AnnotatedIntervalToSegmentVariantContextConverter.convert(feature, referenceContext));

        // Get the correct reference for B37/HG19 compliance:
        // This is necessary because of the variant transformation that gets applied in VariantWalkerBase::apply.
        final ReferenceContext correctReferenceContext = funcotatorEngine.getCorrectReferenceContext(segmentVariantContext, referenceContext);

        // funcotate
        //  The resulting funcotation map should only have one transcript ID (which is the "no transcript" ID).
        final FuncotationMap funcotationMap = funcotatorEngine.createFuncotationMapForSegment(segmentVariantContext, correctReferenceContext, featureContext);

        // write the variant context
        outputRenderer.write(segmentVariantContext, funcotationMap);
    }

    @Override
    public File getDrivingFeatureFile() {
        return segmentFile;
    }

    @Override
    public Object onTraversalSuccess() {
        return true;
    }

    @Override
    public void closeTool() {
        if ( funcotatorEngine != null) {
            funcotatorEngine.close();
        }

        if ( outputRenderer != null ) {
            outputRenderer.close();
        }
    }

    @Override
    protected <T extends Feature> SimpleInterval makeFeatureInterval(T feature) {
        if (funcotatorArgs.referenceVersion.equals("hg19")) {
            return new SimpleInterval(FuncotatorUtils.convertB37ContigToHg19Contig(feature.getContig()), feature.getStart(), feature.getEnd());
        } else {
            return new SimpleInterval(feature);
        }
    }
}
