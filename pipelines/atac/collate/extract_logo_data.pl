# action:
#     parse and aggregate Tn5 cleavage site sequences
# input:
#     19-mer site sequences on STDIN
# outputs:
#     logo plot data on STDOUT

use strict;
use warnings;

use constant {
    KMER_LENGTH     => 19,
    MAX_INDEX       => 18,
    ROUND_FACTOR    => 10000
};

# preassemble list of all possible site sequences, load into hash

# variables
my %baseIndices = (A => 0, C => 1, G => 2, T => 3);
my @baseCounts = (0) x KMER_LENGTH;
@baseCounts = ([@baseCounts], [@baseCounts], [@baseCounts], [@baseCounts]);
my $nSites = 0;

# run the site sequences
while(my $siteSeq = <STDIN>){ 
    chomp $siteSeq; 
    $siteSeq or next;
    my @bases = split("", $siteSeq); 
    foreach my $i(0..MAX_INDEX){
        $baseCounts[$baseIndices{$bases[$i]}][$i] += 1;
    }
    $nSites++;

    # mask unused positions; count site sequences
}

# sort site sequences by count, weight to max count

# print the counts
foreach my $baseI(0..3){
    my $delim = "";
    foreach my $i(0..MAX_INDEX){
        my $freq = int($baseCounts[$baseI][$i] / $nSites * ROUND_FACTOR + 0.5) /  ROUND_FACTOR;
        print "$delim$freq";
        $delim or $delim = "\t";
    }
    print "\n";
}
