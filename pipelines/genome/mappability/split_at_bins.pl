# action:
#     split mappability runs at bin boundaries
# input:
#     sorted BED3 on STDIN of left-aligned kmer start positions representing mappable runs
# outputs:
#     sorted BED3 on STDOUT with:
#       runs split at bin boundaries
#       run ends extended by KMER_LENGTH to support GC assessment in next stream step

# -----*=====|---------- * = kmer start, k = 6, | = bin boundary
# ------*====|=---------
# -------*===|==--------
# --------*==|===-------
# ---------*=|====------
# ----------*|=====-----
# -----------|*=====----
# -----------|-*=====---
# -----------|--*=====--

# -----******|***------- input run of left-aligned kmers start positions

# -----******|=====----- 1st split output run of starts extended to end of rightmost start
# -----------|***=====-- 2nd split output run

# -----******|---------- 1st run of starts after downstream code calculates GC content and unpads the run end position
# -----------|***------- 2nd run

# notice that:
#   bases contributing to runs can cross into the next bin for GC assessment
#   bases near bin boundaries may contribute to two output runs

use strict;
use warnings;

use constant { 
    chrom  => 0,
    start0 => 1,
    end1   => 2,
};

# variables
my $BIN_SIZE    = $ENV{BIN_SIZE};
my $KMER_LENGTH = $ENV{KMER_LENGTH};

# run the mappability runs
while (my $run = <STDIN>) {
    chomp $run;
    my @run = split("\t", $run); 
    my $splitStart0 = $run[start0];
    while ($splitStart0 < $run[end1]) {
        my $binEnd1 = int($splitStart0 / $BIN_SIZE + 1) * $BIN_SIZE;
        my $splitEnd1 = ($run[end1] < $binEnd1) ? $run[end1] : $binEnd1;
        print join("\t", $run[chrom], $splitStart0, $splitEnd1 + $KMER_LENGTH - 1), "\n";
        $splitStart0 = $splitEnd1;
    }
}
