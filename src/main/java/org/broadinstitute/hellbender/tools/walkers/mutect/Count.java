package org.broadinstitute.hellbender.tools.walkers.mutect;

public class Count {
    private final int altCount;
    private final int refCount;

    public Count(final int altCount, final int refCount) {
        this.altCount = altCount;
        this.refCount = refCount;
    }

    public int getAltCount() { return altCount; }

    public int getRefCount() { return refCount; }

    public double getAF() { return (double) altCount / (altCount + refCount); }

    public int getTotal() { return altCount + refCount; }
}
