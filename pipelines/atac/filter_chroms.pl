# action:
#     remove SAM alignments on non-nuclear chroms
# input:
#     headerless SAM stream on STDIN
# outputs:
#     filtered SAM stream on STDOUT

use strict;
use warnings;

# constants
use constant { 
    RNAME => 2
};

# variables
my %allowed;
foreach my $chrom(1..22){
    $allowed{"chr$chrom-$ENV{PRIMARY_GENOME}"}  = 1; # composite chrom names
    $allowed{"chr$chrom-$ENV{SPIKE_IN_GENOME}"} = 1;
}
foreach my $chrom(qw(X Y 2L 2R 3L 3R)){
    $allowed{"chr$chrom-$ENV{PRIMARY_GENOME}"}  = 1;
    $allowed{"chr$chrom-$ENV{SPIKE_IN_GENOME}"} = 1;
}

# run the alignment data
while(my $aln = <STDIN>){ 
    my @f = split("\t", $aln, 4); 
    $allowed{$f[RNAME]} and print join("\t", @f); # skip non-nuclear chroms
}
