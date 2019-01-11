package org.broadinstitute.hellbender.tools.walkers.contamination;


import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary="",
        oneLineSummary = "",
        programGroup = DiagnosticsAndQCProgramGroup.class
)

public class GatherPileupSummaries extends CommandLineProgram {
    // TODO: is it wasteful to ask for a full reference when dictionary will do?
    @Argument(fullName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME, doc = "")
    final File reference = null;

    @Argument(fullName = StandardArgumentDefinitions.INPUT_SHORT_NAME, doc = "")
    final List<File> input = null;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "")
    final File output = null;

    SAMSequenceDictionary sequenceDictionary = null;

    @Override
    protected void onStartup(){
        sequenceDictionary = new ReferenceFileSource(reference.toPath()).getSequenceDictionary();
    }


    @Override
    protected Object doWork() {
        final List<File> nonEmptyFiles = removeEmptyFiles(input);
        Collections.sort(nonEmptyFiles, new PileupSummaryFileComparator(sequenceDictionary));
        PileupSummary.writeToFile(nonEmptyFiles, output);
        return String.format("Successfully merged %d samples", nonEmptyFiles.size());
    }

    private List<File> removeEmptyFiles(final List<File> list){
        final List<File> nonEmptyList = list.stream().filter(f -> PileupSummary.readFromFile(f).getRight().size() > 0)
                .collect(Collectors.toList());
        if (nonEmptyList.size() < list.size()){
            logger.info(String.format("Removed %d empty samples", list.size() - nonEmptyList.size()));
        }
        return nonEmptyList;
    }

    private class PileupSummaryFileComparator implements Comparator<File> {
        final SAMSequenceDictionary sequenceDictionary;
        final List<String> contigsInOrder;

        private PileupSummaryFileComparator(final SAMSequenceDictionary sequenceDictionary){
            this.sequenceDictionary = sequenceDictionary;
            contigsInOrder = sequenceDictionary.getSequences().stream()
                    .map(ssr -> ssr.getSequenceName())
                    .collect(Collectors.toList());
        }
        /**
         * This method makes two assumptions:
         *   1. Each file is already sorted
         *   2. Two files do not overlap (for now)
         */
        @Override
        public int compare(File file1, File file2) {
            final ImmutablePair<String, List<PileupSummary>> pair1 =PileupSummary.readFromFile(file1);
            final ImmutablePair<String, List<PileupSummary>> pair2 =PileupSummary.readFromFile(file2);
            final String sample1 = pair1.getLeft();
            final String sample2 = pair2.getLeft();

            if (! sample1.equals(sample2)){
                throw new UserException(String.format("Cannot merge PileupSummaries from different samples. %s %s"));
            }

            final List<PileupSummary> ps1 = pair1.getRight();
            final PileupSummary first1 = ps1.get(0);

            final List<PileupSummary> ps2 = pair2.getRight();
            final PileupSummary first2 = ps2.get(0);

            // Use Contig Index in case the contig name is e.g. chr2
            final int contigIndex1 = contigsInOrder.indexOf(first1.getContig());
            final int contigIndex2 = contigsInOrder.indexOf(first2.getContig());

            if (contigIndex1 != contigIndex2){
                return contigIndex1 - contigIndex2;
            } else {
                return first1.getStart() - first2.getStart();
            }
        }
    }


}
