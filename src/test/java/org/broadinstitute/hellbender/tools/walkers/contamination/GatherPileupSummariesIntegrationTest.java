package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.walkers.SplitIntervals;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import java.util.List;
import java.util.stream.IntStream;

public class GatherPileupSummariesIntegrationTest extends CommandLineProgramTest {
    @Test
    public void testScatterGather() throws IOException {
        final File intervalDir = createTempDir("intervals");
        final File pileupSummaryDir = createTempDir("pileupSummary");
        final String scatterCount = "20";

        // Step 1: SplitIntervals
        runCommandLine(Arrays.asList(
                "-R", b37_reference_20_21,
                "-L", NA12878_20_21_covered_regions,
                "-O", intervalDir.getAbsolutePath(),
                "-" + SplitIntervals.SCATTER_COUNT_SHORT_NAME, scatterCount),
                SplitIntervals.class.getSimpleName());

        // Step 2: GetPileupSummaries
        final File[] intervals = intervalDir.listFiles();
        IntStream.range(0, intervals.length).forEach(i ->
            runCommandLine(Arrays.asList(
                    "-I", NA12878_20_21_WGS_bam,
                    "-V", thousandGenomes,
                    "-L", intervals[i].getAbsolutePath(),
                    "-O", pileupSummaryDir.getAbsolutePath() + "/" + i + ".tsv",
                    "-" + GetPileupSummaries.MAX_SITE_AF_SHORT_NAME, "0.9"),
                    GetPileupSummaries.class.getSimpleName())
        );

        // Step 3: GatherPileupSummaries
        final File combinedPileupSummary = createTempFile("combined", "tsv");
        final File[] pileupSummaries = pileupSummaryDir.listFiles();
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addArgument(StandardArgumentDefinitions.REFERENCE_SHORT_NAME, b37_reference_20_21);
        args.addArgument(StandardArgumentDefinitions.OUTPUT_SHORT_NAME, combinedPileupSummary.getAbsolutePath());
        Arrays.stream(pileupSummaries).forEach(f -> {
            args.addArgument(StandardArgumentDefinitions.INPUT_SHORT_NAME, f.getAbsolutePath());
        });

        runCommandLine(args.getArgsList(), GatherPileupSummaries.class.getSimpleName());

        final List<PileupSummary> combined = PileupSummary.readFromFile(combinedPileupSummary).getRight();

        Assert.assertTrue(IntStream.range(0, combined.size()-1).allMatch(i ->
                combined.get(i).getStart() <= combined.get(i+1).getStart() ||
                        Integer.valueOf(combined.get(i).getContig()) < Integer.valueOf(combined.get(i+1).getContig())));
    }
}