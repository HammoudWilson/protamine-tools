# action:
#     calculate nuc (nucleosome) and subnuc (subnucleosome) base coverage for all bins on a chromosome
# input:
#     stream of unpacked insert endpoints in BED3 format "chrom\tstart0\tend1"
# outputs:
#     see format below, one line per (moving window of) bin(s)

use strict;
use warnings;

# constants
use constant { 
    CHROM    => 0,
    START0   => 1,
    END1     => 2,
};

# variables
my $BIN_SIZE     = $ENV{AI_BIN_SIZE};
my $MIN_NUC_SIZE = $ENV{MIN_NUC_SIZE};
my $MAX_NUC_SIZE = $ENV{MAX_NUC_SIZE};
my $CHROM        = $ENV{CHROM};
my %baseCounts = map { $_ => [(0) x ($ENV{MAX_BIN_I0} + 100)] } qw(nuc subnuc); # padded on right for moving window
my $binsPerWindow = $MIN_NUC_SIZE / $BIN_SIZE - 1; # will use in moving window below; enforces a 125bp span that must look like a positioned nucleosome

# run the alignment data
while(<STDIN>){ 
    chomp; 
    my @f = split("\t"); 
    my $fragSize = $f[END1] - $f[START0];
    my $type;
    if($fragSize < $MIN_NUC_SIZE){
        $type = 'subnuc';
    } elsif($fragSize <= $MAX_NUC_SIZE){
        $type = 'nuc';
    }
    $type or next; # skip fragments that are more than mononucleosome length
    my $startBin0 = int($f[START0]     / $BIN_SIZE);
    my $endBin0   = int(($f[END1] - 1) / $BIN_SIZE);
    if($startBin0 == $endBin0){
        $baseCounts{$type}[$startBin0] += $fragSize;
    } else {
        my $startBinOffset = $f[START0]     % $BIN_SIZE;
        my $endBinOffset   = ($f[END1] - 1) % $BIN_SIZE;
        $baseCounts{$type}[$startBin0] += $BIN_SIZE - $startBinOffset;
        $baseCounts{$type}[$endBin0]   += $endBinOffset;
        if($endBin0 > $startBin0 + 1){
            foreach my $bin0 (($startBin0 + 1)..($endBin0 - 1)){
                $baseCounts{$type}[$bin0] += $BIN_SIZE;
            }
        }
    }
}

sub sum{
    my $sum = 0;
    foreach (@_) { $sum += $_ }
    return $sum;
}

# print the counts
print(join("\t", 
    "chrom",
    "start0",
    "nuc_bin", 
    "subnuc_bin",
    "nuc_mean",
    "subnuc_mean",
    "bias"
), "\n");
foreach my $binI0(0..$ENV{MAX_BIN_I0}){
    my @i = $binI0..($binI0 + $binsPerWindow - 1);
    my $nucSum    = sum(@{$baseCounts{nuc   }}[@i]);
    my $subnucSum = sum(@{$baseCounts{subnuc}}[@i]);
    my $sum = $nucSum + $subnucSum;
    print join("\t", 
        $CHROM,
        $binI0 * $BIN_SIZE,
        $baseCounts{nuc   }[$binI0],
        $baseCounts{subnuc}[$binI0],
        $nucSum    / $binsPerWindow,
        $subnucSum / $binsPerWindow,
        $sum > 0 ? ($nucSum - $subnucSum) / $sum : 0
    ), "\n"
}
