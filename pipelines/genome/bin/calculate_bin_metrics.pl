# action:
#     calculate bin metrics from mappability intersection
#       fraction of mappable kmers in each bin
#       weighted average of the GC content of those mappable kmers
#     thus, bin mappability and GC content reflects just the mappable inserts and may differ by kmer
# input:
#     bedtools intersect of BED3 bins to mappability runs previously split at bin boundaries
# outputs:
#     grouped and aggregated bin-level metrics as per above

use strict;
use warnings;

use constant { 
    chrom_bin   => 0, # typically 1kb bins, might have multiple runs per bin
    start0_bin  => 1,
    end1_bin    => 2,
    chrom_run   => 3, # one run per bin per line ...
    start0_run  => 4,
    end1_run    => 5,
    run_frac_gc => 6, # ... with GC metrics per run
    n_bases_overlap => 7, # ... and the number of bases overlap between bin and mappability run
};

# variables
my $BIN_SIZE = $ENV{BIN_SIZE};
my ($prevBinKey, @intersects);

# thread through the bin-to-run intersections
while (my $intersect = <STDIN>) {
    chomp $intersect;
    my @intersect = split("\t", $intersect); 
    my $binKey = join("\t", @intersect[chrom_bin, start0_bin]);
    if($prevBinKey and $prevBinKey ne $binKey) {
        processBin();
        @intersects = ();
    }
    $prevBinKey = $binKey;
    push @intersects, \@intersect;
}
processBin();

# calculate bin metrics
# print just the fraction of mappable kmers in each bin and the weighted average of the GC content
sub processBin {
    if($intersects[0][chrom_run] ne '.') { # bedtools intersect -wao identifies no b spans as chr == "."
        my $nMappablePos  = 0;
        my $sumWeightedGC = 0;
        foreach my $intersect (@intersects) {
            $nMappablePos += $$intersect[n_bases_overlap];
            $sumWeightedGC += $$intersect[run_frac_gc] * $$intersect[n_bases_overlap];
        }
        print join("\t", 
            sprintf("%.4f", $nMappablePos  / $BIN_SIZE), 
            sprintf("%.4f", $sumWeightedGC / $nMappablePos),
        ), "\n";
    } else {
        print join("\t", 0, 0), "\n";
    }
}
