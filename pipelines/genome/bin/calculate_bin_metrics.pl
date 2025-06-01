# action:
#     calculate bin metrics from mappability intersection at a given insert size level
#       fraction of mappable left-aligned insert positions in each bin
#       weighted average of the GC content of those mappable inserts
#     thus, values reflect just the mappable inserts in a bin and may differ by insert size level
# input:
#     bedtools intersect of BED3 bins to mappability runs previously split at bin boundaries
# outputs:
#     grouped and aggregated bin-level metrics as mpp_bin_N\tgc_bin_N

use strict;
use warnings;

use constant { 
    CHROM_BIN   => 0, # typically 1kb bins, might have multiple mappability runs per bin
    START0_BIN  => 1,
    END1_BIN    => 2,
    CHROM_RUN   => 3, # one mappability run per bin per line ...
    START0_RUN  => 4,
    END0_RUN    => 5,
    RUN_FRAC_GC => 6, # ... with GC metrics per run
    N_BASES_OVERLAP => 7, # ... and the number of bases overlap between bin and mappability run
};

# variables
my $BIN_SIZE = $ENV{BIN_SIZE};
my ($prevBinKey, @intersects);

# thread through the bin-to-run intersections
while (my $intersect = <STDIN>) {
    chomp $intersect;
    my @intersect = split("\t", $intersect); 
    my $binKey = join("\t", @intersect[CHROM_BIN, START0_BIN]);
    if($prevBinKey and $prevBinKey ne $binKey) {
        processBin();
        @intersects = ();
    }
    $prevBinKey = $binKey;
    push @intersects, \@intersect;
}
processBin();

# calculate bin metrics
#   fraction of mappable start positions in each bin
#   weighted average of the GC content
sub processBin {
    if($intersects[0][CHROM_RUN] eq '.') { # bedtools intersect -wao identifies no b spans as chr == "."
        print join("\t", 0, 0), "\n";
    } else {
        my $nMappablePos  = 0;
        my $sumWeightedGC = 0;
        foreach my $intersect (@intersects) {
            $nMappablePos += $$intersect[N_BASES_OVERLAP];
            $sumWeightedGC += $$intersect[RUN_FRAC_GC] * $$intersect[N_BASES_OVERLAP];
            # remember, mappability runs are split at bin boundaries, so this weighting is correct
            # see genome/mappability for details
        }
        print join("\t", 
            sprintf("%.4f", $nMappablePos  / $BIN_SIZE), 
            sprintf("%.4f", $sumWeightedGC / $nMappablePos)
        ), "\n";
    }
}
