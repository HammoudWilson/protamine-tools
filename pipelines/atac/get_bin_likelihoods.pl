# action:
#     unpack the distinct insert end1s for every insert start1
#     calculate the log likelihood of each insert for histone- and protamine-associated models
#     sum the log likelihoods for each bin over all inserts overlapping it
# input:
#     stream of distinct insert endpoints in format "start1\tend1_1[,end1_2,...]"
# outputs:
#     LL_histone, LL_protamine, nInserts for each bin on STDOUT
#     to support downstream parsing for bin-level NRLL calculation or HMM

use strict;
use warnings;

# constants
use constant { 
    START1 => 0,
    END1S  => 1
};

# variables
my $BIN_SIZE = $ENV{BIN_SIZE};
my @nInserts = (0) x ($ENV{MAX_BIN_I0} + 1);
my %logLikelihoods = (
    histone   => [@nInserts],
    protamine => [@nInserts]
);

# load emission probabilities
my %eps;
open(my $fh, $ENV{EMISS_PROBS_FILE});
while(<$fh>){
    chomp;
    my @f = split("\t");
    push @{$eps{histone}},   $f[0];
    push @{$eps{protamine}}, $f[1];
}
close($fh);

# run the alignment data
while(<STDIN>){ 
    chomp; 
    my @f = split("\t"); 
    foreach my $end1(split(",", $f[END1S])){ 
        my $insertSize0 = $end1 - $f[START1];
        my $LL_histone   = $eps{histone}[  $insertSize0];
        my $LL_protamine = $eps{protamine}[$insertSize0];
        my $startBin0   = int(($f[START1] - 1) / $BIN_SIZE);
        my $endBin0     = int(($end1      - 1) / $BIN_SIZE);
        $logLikelihoods{histone  }[$startBin0] += $LL_histone;
        $logLikelihoods{protamine}[$startBin0] += $LL_protamine;
        $logLikelihoods{histone  }[$endBin0]   += $LL_histone;
        $logLikelihoods{protamine}[$endBin0]   += $LL_protamine;
        $nInserts[$startBin0]++; # OK to double count both LL and nInserts for single-bin inserts since normalizing below
        $nInserts[$endBin0]++;
    }
}

# print the counts
foreach my $binI0(0..$ENV{MAX_BIN_I0}){
    print join("\t", 
        $logLikelihoods{histone  }[$binI0],
        $logLikelihoods{protamine}[$binI0],
        $nInserts[$binI0]
    ), "\n";
}
