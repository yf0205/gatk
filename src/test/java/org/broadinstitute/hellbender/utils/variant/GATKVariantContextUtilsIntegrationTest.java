package org.broadinstitute.hellbender.utils.variant;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.util.TabixUtils;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderVersion;
import java.io.File;
import java.io.IOException;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public final class GATKVariantContextUtilsIntegrationTest extends GATKBaseTest {

    // Same test as in GATKVariantContextUtilsUnitTest, but writing to a real bucket.
    @Test(groups = {"bucket"})
    public void testCreateVcfWriterOnNio() throws IOException {
        final File inputGZIPFile = new File(
            publicTestDir + "org/broadinstitute/hellbender/engine/8_mutect2_sorted.vcf.gz");
        final String dest = BucketUtils.getTempFilePath(
            getGCPTestStaging() +"testCreateVcfWriterOnNio", ".vcf.gz");
        final Path outputGZIP = BucketUtils.getPathOnGcs(dest);
        final Path tabixIndex = outputGZIP.resolveSibling(outputGZIP.getFileName().toString() + TabixUtils.STANDARD_INDEX_EXTENSION);
        long recordCount = 0;

        try (final VariantContextWriter vcfWriter = GATKVariantContextUtils.createVCFWriter(
            outputGZIP,
            null,
            false,
            Options.INDEX_ON_THE_FLY);
            final FeatureReader<VariantContext> inputFileReader =
                AbstractFeatureReader
                    .getFeatureReader(inputGZIPFile.getAbsolutePath(), new VCFCodec(), false)) {
            vcfWriter.writeHeader((VCFHeader) inputFileReader.getHeader());
            final Iterator<VariantContext> it = inputFileReader.iterator();
            while (it.hasNext()) {
                vcfWriter.add(it.next());
                recordCount++;
            }
        }

        // make sure we got an output and tabix index
        Assert.assertTrue(Files.exists(outputGZIP));
        Assert.assertTrue(Files.exists(tabixIndex));
        Assert.assertTrue(Files.size(tabixIndex) > 0);

        // make sure we got an output and queryable index
        long roundTripRecordCount = 0;
        try (final VCFFileReader outputFileReader = new VCFFileReader(outputGZIP, true)) {
            final Iterator<VariantContext> it = outputFileReader.query(new SimpleInterval("chr6", 1, 999999999));
            while (it.hasNext()) {
                it.next();
                roundTripRecordCount++;
            }
        }
        Assert.assertEquals(roundTripRecordCount, recordCount);
    }
}
