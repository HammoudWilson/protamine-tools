# action:
#     unpack the distinct insert end1s for every insert start1
#     count the number of insert endpoints in each bin, establishing counts for:
#         all inserts of any size
#         intermediate-size inserts >= 65 bp and <= 125 bp, typical of bins during spermatid elongation and protamine replacement
# input:
#     stream of distinct insert endpoints in format "start1\tend1_1[,end1_2,...]"
# outputs:
#     fractional independent read pair counts (in increments of 0.5) for each possible bin on chromosome
#     a read that crosses a bin boundary will be counted as 0.5 in each bin
#     due to filtering upstream, no insert will span more than 2 bins
#     at present the method does not take into account the fraction of a read in each bin
#         (unlikely to be important since bins are large relative to inserts)

use strict;
use warnings;

# constants
use constant { 
    ALL_INSERTS => "all_inserts",
    MIN_INTERMEDIATE_SIZE_BP => 65,
    MAX_INTERMEDIATE_SIZE_BP => 125,
    START1 => 0,
    END1S  => 1
};

# variables
my $INSERT_TYPE = $ENV{INSERT_TYPE};
my $BIN_SIZE = $ENV{BIN_SIZE};
my @counts = (0) x ($ENV{MAX_BIN_I0} + 1);

# run the alignment data
while(<STDIN>){ 
    chomp; 
    my @f = split("\t"); 
    foreach my $end1(split(",", $f[END1S])){ 
        my $insertSize = $end1 - $f[START1] + 1;
        if(
            $INSERT_TYPE eq ALL_INSERTS or 
           (
                $insertSize >= MIN_INTERMEDIATE_SIZE_BP and 
                $insertSize <= MAX_INTERMEDIATE_SIZE_BP
            )
        ){
                $counts[int(($f[START1] - 1) / $BIN_SIZE)] += 0.5; # 0.5 for each end, thus, 1.0 for each read pair
                $counts[int(($end1      - 1) / $BIN_SIZE)] += 0.5;
        }
    }
}

# print the counts
foreach my $binI0(0..$ENV{MAX_BIN_I0}){
    print $counts[$binI0], "\n";
}
