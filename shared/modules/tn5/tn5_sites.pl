# support for assessing Tn5 cleavage sites

use strict;
use warnings;

# variables
our @MAPPABILITY_SIZE_LEVELS = split(/\s+/, $ENV{MAPPABILITY_SIZE_LEVELS});
our @TN5_PREFERENCE_POSITIONS = split(/\s+/, $ENV{TN5_PREFERENCE_POSITIONS});
our $N_TN5_PREFERENCE_POSITIONS = $ENV{N_TN5_PREFERENCE_POSITIONS};
our $nTn5Sites = 4 ** $N_TN5_PREFERENCE_POSITIONS;
our $PRIMARY_GENOME  = $ENV{PRIMARY_GENOME};
our $SPIKE_IN_GENOME = $ENV{SPIKE_IN_GENOME};

# load the possible Tn5 cleavage site nonamers
our (@tn5Sites, %tn5SiteIs, %tn5SiteCounts);
open my $inH, "<", $ENV{TN5_KMERS_FILE} or die "collate_count_inserts error: cannot open $ENV{TN5_KMERS_FILE} for reading\n";
my $siteI = 0;
while(<$inH>){
    chomp;
    push @tn5Sites, $_;
    $tn5SiteIs{$_} = $siteI;
    $siteI++;
}
close $inH;

# initialize Tn5 nonamer counting by insert size level
our (@insertSizeLevelIs);
my $islI = 0;
foreach my $insertSizeLevel(@MAPPABILITY_SIZE_LEVELS){
    $insertSizeLevelIs[$insertSizeLevel] = $islI;
    $tn5SiteCounts{$PRIMARY_GENOME }[$islI] = [(0) x $nTn5Sites];
    $tn5SiteCounts{$SPIKE_IN_GENOME}[$islI] = [(0) x $nTn5Sites];
    $islI++;
}

# increment the count of a Tn5 cleavage site, return the parsed discontiguous nonamer
sub countTn5Site {
    my ($genome, $islI, $ucSeq13) = @_;

    # extract just TN5_PREFERENCE_POSITIONS from $seq13
    # thus, the input 13-mer is reduced a discontiguous 9-mer
    # note that nonamers are always on the top reference strand at either end of an insert
    # caller is responsible for
    #   providing ucSeq13 in the proper register depending on the cleavage position and strand
    #   ensuring that ucSeq13 is uppercase
    #    0  v
    # ---**|*--***--* **---
    # ---** *--***--*|**---
    #               ^  0
    my $nonamer = '';
    foreach my $pos (@TN5_PREFERENCE_POSITIONS) {
        $nonamer .= substr($ucSeq13, $pos, 1);
    }

    # drop any nonamers that contain N bases, which will occur when assessing genomic sites
    exists $tn5SiteIs{$nonamer} or return;

    # count expected nonamers
    $tn5SiteCounts{$genome}[$islI][$tn5SiteIs{$nonamer}]++;

    # return nonamer for recording with insert spans
    $nonamer; 
}

# print the Tn5 cleavage site nonamer counts to STDOUT for reading in R
# site sequence not printed, guaranteed to be in the same order as $TN5_KMERS_FILE
sub printTn5SiteCounts {
    print join("\t", 
        "genome", 
        map { "tn5_" . $_ } @MAPPABILITY_SIZE_LEVELS
    ), "\n";
    foreach my $genome ($PRIMARY_GENOME, $SPIKE_IN_GENOME) {
        foreach my $nonamer(@tn5Sites){
            print join("\t", 
                $genome,
                map { $tn5SiteCounts{$genome}[$_][$tn5SiteIs{$nonamer}] } 0..$#MAPPABILITY_SIZE_LEVELS
            ), "\n";
        }
    }
}

1;
