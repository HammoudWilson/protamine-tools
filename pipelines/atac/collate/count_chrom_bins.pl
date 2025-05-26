# action:
#     unpack the distinct insert end1s for every insert start1
#     count the number of left-justified insert endpoints in each bin
# input:
#     stream of distinct insert endpoints in format "chrom\tstart0\tend1_1[,end1_2,...]"
# outputs:
#     independent read pair counts for each possible bin on chromosome
#     a read that crosses a bin boundary will be counted as 1 in the leftmost bin
#       this choice is made to maximize the accuracy of mappability, as genmap values are left justified
#     due to filtering upstream, no insert will span more than 2 bins
#     expect ~90% of reads to be restricted to a single bin based on bin size and insert size distributions

use strict;
use warnings;

# constants
use constant { 
    CHROM  => 0,
    START0 => 1,
    END1S  => 2
};

# variables
my $BIN_SIZE = $ENV{BIN_SIZE};
my @counts = (0) x ($ENV{MAX_BIN_I0} + 1);

# run the alignment data
while(<STDIN>){ 
    chomp; 
    my @f = split("\t"); 
    foreach my $end1(split(",", $f[END1S])){ 
        $counts[int(($f[START0]) / $BIN_SIZE)]++;
        # $counts[int(($f[START0]) / $BIN_SIZE)] += 0.5; # 0.5 for each end, thus, 1.0 for each read pair
        # $counts[int(($end1- 1)   / $BIN_SIZE)] += 0.5;
    }
}

# print the counts
foreach my $binI0(0..$ENV{MAX_BIN_I0}){
    print $counts[$binI0], "\n";
}
