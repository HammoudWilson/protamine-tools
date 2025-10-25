# action:
#     calculate the log likelihood of each insert for histone- and protamine-associated models
#     sum the log likelihoods for each bin over all inserts mapping to it
# input:
#     stream of distinct inserts in format chrom,start0,end1,insertSizeLevel
# outputs:
#     LL_histone, LL_protamine, nInserts for each bin on STDOUT
#     to support downstream parsing for bin-level NRLL calculation or HMM

use strict;
use warnings;

# constants
use constant { 
    CHROM             => 0,
    START0            => 1,
    END1              => 2,
    INSERT_SIZE_LEVEL => 3,
};

# variables
my $BIN_SIZE = $ENV{BIN_SIZE};
my @MAPPABILITY_SIZE_LEVELS = split(/\s+/, $ENV{MAPPABILITY_SIZE_LEVELS});
my (%logLikelihoods, %nInserts);

# load emission probabilities
my %eps_state_isl;
open(my $fh, $ENV{EMISS_PROBS_FILE});
while(<$fh>){
    chomp;
    my @f = split("\t");
    push @{$eps_state_isl{histone}},   $f[0];
    push @{$eps_state_isl{protamine}}, $f[1];
}
close($fh);

# initialize insert size level indices
my (@insertSizeLevelIs);
my $islI = 0;
foreach my $insertSizeLevel(@MAPPABILITY_SIZE_LEVELS){
    $insertSizeLevelIs[$insertSizeLevel] = $islI;
    $islI++;
}

# run the individual inserts
while(<STDIN>){ 
    chomp; 
    my @ins = split("\t"); 
    my $islI = $insertSizeLevelIs[$ins[INSERT_SIZE_LEVEL]];
    my $LL_histone   = $eps_state_isl{histone}[  $islI];
    my $LL_protamine = $eps_state_isl{protamine}[$islI];
    my $startBin0 = int($ins[START0] / $BIN_SIZE) * $BIN_SIZE;
    my $binKey = "$ins[CHROM]:$startBin0";
    $logLikelihoods{histone  }{$binKey} += $LL_histone;
    $logLikelihoods{protamine}{$binKey} += $LL_protamine;
    $nInserts{$binKey}++;
}

# print the LL and total counts for all bins; use zero for bins without inserts
open my $binH, "-|", "zcat $ENV{GENOME_BINS_BED}" 
    or die "get_bin_likelihoods error: cannot open $ENV{GENOME_BINS_BED} for reading\n";
my $header = <$binH>;
while(<$binH>){
    chomp;
    my @bin = split("\t");
    my $binKey = "$bin[CHROM]:$bin[START0]";
    print join("\t", 
        $logLikelihoods{histone  }{$binKey} || 0,
        $logLikelihoods{protamine}{$binKey} || 0,
        $nInserts{$binKey} || 0
    ), "\n";
}
close $binH;
