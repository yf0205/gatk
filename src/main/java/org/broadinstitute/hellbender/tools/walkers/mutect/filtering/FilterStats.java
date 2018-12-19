package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.Collection;

public class FilterStats {
    private final String filterName;
    private final double threshold;
    private final double expectedNumFPs;
    private final int numPassingVariants;
    private final double expectedFPR;

    public FilterStats(final String filterName, final double threshold, final double expectedNumFPs,
                       final int numPassingVariants, final double expectedFPR){
        this.filterName = filterName;
        this.threshold = threshold;
        this.expectedNumFPs = expectedNumFPs;
        this.numPassingVariants = numPassingVariants;
        this.expectedFPR = expectedFPR;
    }

    public String getFilterName() { return filterName; }

    public double getExpectedNumFPs() { return expectedNumFPs; }

    public int getNumPassingVariants() { return numPassingVariants; }

    public double getThreshold() { return threshold; }

    public double getExpectedFPR() { return expectedFPR; }

    private enum M2FilterStatsTableColumn {
        FILTER_NAME("filter_name"),
        THRESHOLD("threshold"),
        EXPECTED_FALSE_POSITIVES("expected_fps"),
        EXPECTED_FALSE_POSITIVE_RATE("expected_fpr"),
        NUM_PASSING_VARIANTS("num_passing_variants");

        private String columnName;

        M2FilterStatsTableColumn(final String columnName) {
            this.columnName = columnName;
        }

        @Override
        public String toString() { return columnName; }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static class Mutect2FilterStatsWriter extends TableWriter<FilterStats> {
        private Mutect2FilterStatsWriter(final File output) throws IOException {
            super(output, M2FilterStatsTableColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final FilterStats stats, final DataLine dataLine) {
            dataLine.set(M2FilterStatsTableColumn.FILTER_NAME.toString(), stats.getFilterName())
                    .set(M2FilterStatsTableColumn.THRESHOLD.toString(), stats.getThreshold())
                    .set(M2FilterStatsTableColumn.EXPECTED_FALSE_POSITIVES.toString(), stats.getExpectedNumFPs())
                    .set(M2FilterStatsTableColumn.EXPECTED_FALSE_POSITIVE_RATE.toString(), stats.getExpectedFPR())
                    .set(M2FilterStatsTableColumn.NUM_PASSING_VARIANTS.toString(), stats.getNumPassingVariants());
        }
    }

    public static void writeM2FilterSummary(final Collection<FilterStats> filterStats, final File outputTable) {
        try (Mutect2FilterStatsWriter writer = new Mutect2FilterStatsWriter(outputTable)) {
            writer.writeAllRecords(filterStats);
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", outputTable), e);
        }
    }

}
