package org.broadinstitute.hellbender.tools.funcotator;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalCollection;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class FuncotateSegmentsIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_SUB_DIR = toolsTestDir + "/funcotator/";
    private static final String SIMPLE_TEST_FILE = TEST_SUB_DIR + "simple.seg";
    private static final String SIMPLE_TEST_CNTN4_FILE = TEST_SUB_DIR + "simple_cntn4_overlap.seg";
    private static final String REF = b37Reference;
    private static final String DS_PIK3CA_DIR  = largeFileTestDir + "funcotator" + File.separator + "small_ds_pik3ca" + File.separator;
    // This has transcripts with multiple gene names...
    private static final String DS_CNTN4_DIR  = toolsTestDir + "funcotator" + File.separator + "small_cntn4_ds" + File.separator;

    @Test
    public void testSimpleNoOverlap() throws IOException {
        final File outputFile = File.createTempFile("funcotatesegs_simple", ".seg");

        final ArgumentsBuilder arguments = new ArgumentsBuilder();
        arguments.add("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.add(SIMPLE_TEST_FILE);
        arguments.add("--" + FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME);
        arguments.add(FuncotatorArgumentDefinitions.OutputFormatType.SEG);
        arguments.add("--" + StandardArgumentDefinitions.REFERENCE_LONG_NAME);
        arguments.add(REF);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        arguments.add("--" + FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME);
        arguments.add("hg19");
        arguments.addArgument(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, DS_PIK3CA_DIR);

        runCommandLine(arguments);

        final AnnotatedIntervalCollection collection = AnnotatedIntervalCollection.create(outputFile.toPath(), null);
        Assert.assertEquals(collection.getRecords().size(), 3);
        Assert.assertTrue(collection.getRecords().stream().allMatch(r -> r.hasAnnotation("genes")));
        Assert.assertTrue(collection.getRecords().stream().allMatch(r -> r.getAnnotationValue("genes").equals("")));
    }

    /**
     * Very dependent on the data in "simple_cntn4_overlap.seg"
     */
    @DataProvider
    public Object[][] cntn4GroundTruth() {
        return new Object[][] {
                {
                            Arrays.asList(
                                    "CNTN4, CNTN4-AS2", "CNTN4", ""
                            )
                }
        };
    }

    @Test(dataProvider = "cntn4GroundTruth")
    public void testSimpleMultipleGenesOverlap(List<String> gtGenesValues) throws IOException {
        final File outputFile = File.createTempFile("funcotatesegs_simple_cntn4", ".seg");

        final ArgumentsBuilder arguments = new ArgumentsBuilder();
        arguments.add("--" + CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME);
        arguments.add(SIMPLE_TEST_CNTN4_FILE);
        arguments.add("--" + FuncotatorArgumentDefinitions.OUTPUT_FORMAT_LONG_NAME);
        arguments.add(FuncotatorArgumentDefinitions.OutputFormatType.SEG);
        arguments.add("--" + StandardArgumentDefinitions.REFERENCE_LONG_NAME);
        arguments.add(REF);
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        arguments.add("--" + FuncotatorArgumentDefinitions.REFERENCE_VERSION_LONG_NAME);
        arguments.add("hg19");
        arguments.addArgument(FuncotatorArgumentDefinitions.DATA_SOURCES_PATH_LONG_NAME, DS_CNTN4_DIR);

        runCommandLine(arguments);

        final AnnotatedIntervalCollection collection = AnnotatedIntervalCollection.create(outputFile.toPath(), null);
        Assert.assertEquals(collection.getRecords().size(), 3);

        final List<String> testGenesValues = collection.getRecords().stream().map(r -> r.getAnnotationValue("genes")).collect(Collectors.toList());
        Assert.assertEquals(testGenesValues, gtGenesValues);
    }
}
