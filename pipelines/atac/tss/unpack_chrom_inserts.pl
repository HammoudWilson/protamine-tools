# action:
#     unpack the distinct insert end1s for every insert start1
# input:
#     stream of distinct insert endpoints in format "start1\tend1_1[,end1_2,...]"
# outputs:
#     stream of unpacked insert endpoints in BED3 format "chrom\tstart0\tend1"
#     sort is maintained the same as the sorted BAM input

use strict;
use warnings;

# constants
use constant { 
    START1 => 0,
    END1S  => 1
};

# variables
my $CHROM = $ENV{CHROM};

# run the alignment data
while(<STDIN>){ 
    chomp; 
    my @f = split("\t"); 
    foreach my $end1(split(",", $f[END1S])){ 
        print join("\t", $CHROM, $f[START1] - 1, $end1)."\n";
    }
}
