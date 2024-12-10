# action:
#     unpack the distinct insert end1s for every insert start1
#     calculate the log likelihood of each insert for nucleosome and subnucleosome models
#     sum the log likelihoods for each bin over all inserts overlapping it
# input:
#     stream of distinct insert endpoints in format "start1\tend1_1[,end1_2,...]"
# outputs:
#     LL_nuc, LL_subnuc, nInserts for each bin on STDOUT
#     formatted to support segment.pl and calculating NRLL for each bin and called HMM segment

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
    nucleosomal    => [@nInserts],
    subnucleosomal => [@nInserts]
);

# load emission probabilities
my %eps;
open(my $fh, $ENV{EMISS_PROBS_FILE});
while(<$fh>){
    chomp;
    my @f = split("\t");
    push @{$eps{nucleosomal}},    $f[0];
    push @{$eps{subnucleosomal}}, $f[1];
}
close($fh);

# run the alignment data
while(<STDIN>){ 
    chomp; 
    my @f = split("\t"); 
    foreach my $end1(split(",", $f[END1S])){ 
        my $insertSize0 = $end1 - $f[START1];
        my $LL_nuc      = $eps{nucleosomal}[   $insertSize0];
        my $LL_subnuc   = $eps{subnucleosomal}[$insertSize0];
        my $startBin0   = int(($f[START1] - 1) / $BIN_SIZE);
        my $endBin0     = int(($end1      - 1) / $BIN_SIZE);
        $logLikelihoods{nucleosomal   }[$startBin0] += $LL_nuc;
        $logLikelihoods{subnucleosomal}[$startBin0] += $LL_subnuc;
        $logLikelihoods{nucleosomal   }[$endBin0]   += $LL_nuc;
        $logLikelihoods{subnucleosomal}[$endBin0]   += $LL_subnuc;
        $nInserts[$startBin0]++; # OK to double count both LL and nInserts for single-bin inserts since normalizing below
        $nInserts[$endBin0]++;
    }
}

# print the counts
foreach my $binI0(0..$ENV{MAX_BIN_I0}){
    my $LL_nuc    = $logLikelihoods{nucleosomal   }[$binI0]; # TODO: divide by 2?
    my $LL_subnuc = $logLikelihoods{subnucleosomal}[$binI0];
    my $nInserts  = $nInserts[$binI0];
    print join("\t", 
        join(",", $LL_nuc, $LL_subnuc), # will be replaced with 0|1 HMM state by segment.pl
        $LL_nuc,
        $LL_subnuc,
        $nInserts
        # ,
        # $nInserts > 0 ? ($LL_subnuc - $LL_nuc) / $nInserts : 0
    ), "\n";
}
