# action:
#     process insersection of ATAC-seq reads with TSSs
#     for each TSS that matched overlapped an insert
#         calculate the strand-oriented distance from the TSS to the insert endpoints
# input:
#     stream of from bedtools intersect -wo -a inserts.bed -b tss.bed
# outputs:
#     potentially redundant BED file of inserts with the following columns:
#         chrom, start0, end1, tss_i1, start_to_tss, end_to_tss

use strict;
use warnings;

# constants
use constant { 
    INSERT_CHROM    => 0,
    INSERT_START0   => 1,
    INSERT_END1     => 2,
    TSS_CHROM       => 3,
    TSS_START1      => 4,
    TSS_END1        => 5,
    TSS_STRAND      => 6,
    TSS_I1          => 7, # row number in TSS file
    #--------------------
    TSS_FLANK_LEN   => 1000,
};

# run the alignment data
print(join("\t", 
    "chrom", "start0", "end1", "tss_i1", "start_to_tss", "end_to_tss"
), "\n");
while(<STDIN>){ 
    chomp; 
    my @f = split("\t"); 
    my $tss1 = $f[TSS_START1] + TSS_FLANK_LEN;
    print(join("\t", 
        @f[INSERT_CHROM..INSERT_END1],
        $f[TSS_I1],
        $f[TSS_STRAND] eq "+" ? 
            $f[INSERT_START0] - $tss1 : 
            $tss1 - $f[INSERT_START0],
        $f[TSS_STRAND] eq "+" ? 
            $f[INSERT_END1] - $tss1 : 
            $tss1 - $f[INSERT_END1]
    ), "\n");
}
