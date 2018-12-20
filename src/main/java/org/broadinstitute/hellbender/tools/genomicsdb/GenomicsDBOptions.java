package org.broadinstitute.hellbender.tools.genomicsdb;

import java.nio.file.Path;

/**
 * : Path to a reference. May be null. Needed only for reading from GenomicsDB.
 *      * @param doGnarlyGenotyping indicates whether the GenomicsDB export configuration should be modified for the GnarlyGenotyper
 */
public class GenomicsDBOptions {
    private Path reference;
    private boolean callGenotypes;

    public GenomicsDBOptions() {
        this(null, false);
    }

    public GenomicsDBOptions(final Path reference) {
        this(reference, false);
    }

    public GenomicsDBOptions(final Path reference, final boolean callGenotypes) {
        this.reference = reference;
        this.callGenotypes = callGenotypes;
    }

    public Path getReference() {
        return reference;
    }

    public boolean doCallGenotypes() {
        return callGenotypes;
    }
}
