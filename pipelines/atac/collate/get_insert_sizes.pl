# action:
#     unpack the distinct insert end1s for every insert start1
#     calculate insert size distributions
# input:
#     stream of distinct insert endpoints in format "start1\tend1_1[,end1_2,...]"
# outputs:
#     insert size distribution on STDOUT

use strict;
use warnings;

# constants
use constant { 
    START1 => 0,
    END1S  => 1
};

# variables
my $BIN_SIZE = $ENV{BIN_SIZE};
my @counts = (0) x $BIN_SIZE; # includes insert sizes 1 to $MIN_INSERT_SIZE, but those counts will always be 0

# run the alignment data
while(<STDIN>){ 
    chomp; 
    my @f = split("\t"); 
    foreach my $end1(split(",", $f[END1S])){ 
        $counts[$end1 - $f[START1]]++; # 0-indexed insert size
    }
}

# print the counts
foreach my $insertSize0(0..($BIN_SIZE - 1)){
    print $counts[$insertSize0], "\n";
}
