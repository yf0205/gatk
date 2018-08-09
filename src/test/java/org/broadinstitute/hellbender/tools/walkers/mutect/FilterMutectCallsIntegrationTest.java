package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.annotations.Test;

import java.io.File;

import static org.testng.Assert.*;

@Test(groups = {"variantcalling"})
public class FilterMutectCallsIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testNormalNormal() {
        Utils.resetRandomGenerator();

        final File inputVcf = new File("/Users/davidben/Desktop/dream/4/", "synthetic.challenge.set4.tumor-filtered.vcf");
        //final File inputVcf = new File("/Users/davidben/Desktop/normal-normal", "nn4.vcf");
        final File outputVcf = createTempFile("filtered", ".vcf");
        final String[] createPonArgs = {
                "-V", inputVcf.getAbsolutePath(),
                "-O", outputVcf.getAbsolutePath()
        };

        runCommandLine(createPonArgs);

        int k = 10;

    }

}