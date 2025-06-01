# action:
#     split mappability runs at bin boundaries
# input:
#     sorted BED3 on STDIN of left-aligned insert start positions representing mappable runs
# outputs:
#     sorted BED3 on STDOUT with:
#       runs split at bin boundaries
#       run ends extended by INSERT_SIZE_LEVEL to support GC assessment in next stream step

# -----*=====|---------- * = insert start, insert size level = 6, | = bin boundary
# ------*====|=---------
# -------*===|==--------
# --------*==|===-------
# ---------*=|====------
# ----------*|=====-----
# -----------|*=====----
# -----------|-*=====---
# -----------|--*=====--

# -----******|***------- input run of left-aligned insert start positions

# -----******|=====----- 1st split output run of starts extended to the insert end of the rightmost start
# -----------|***=====-- 2nd split output run

# -----******|---------- 1st run of starts after downstream code calculates GC content and unpads the run end position
# -----------|***------- 2nd run

# ---oo******|***=====oo the same span after padding and masking in preparation for Tn5 site analysis (handled later)

# notice that:
#   bases contributing to runs can cross into the next bin for GC assessment
#   bases near bin boundaries may contribute to two output runs

use strict;
use warnings;

use constant { 
    CHROM  => 0,
    START0 => 1,
    END1   => 2,
};

# variables
my $BIN_SIZE = $ENV{BIN_SIZE};
my $INSERT_SIZE_LEVEL = $ENV{INSERT_SIZE_LEVEL};

# run the mappability runs
while (my $run = <STDIN>) {
    chomp $run;
    my @run = split("\t", $run); 
    my $splitStart0 = $run[START0];
    while ($splitStart0 < $run[END1]) {
        my $binEnd1 = int($splitStart0 / $BIN_SIZE + 1) * $BIN_SIZE;
        my $splitEnd1 = ($run[END1] < $binEnd1) ? $run[END1] : $binEnd1;
        print join("\t", $run[CHROM], $splitStart0, $splitEnd1 + $INSERT_SIZE_LEVEL - 1), "\n";
        $splitStart0 = $splitEnd1;
    }
}
