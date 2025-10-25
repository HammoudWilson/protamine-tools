# action:
#     calculate insert size distributions for primary genome, stratified by insert percent GC
#     these distributions are used to:
#       assess relative insert sizes
#       explore the correlation of insert size and GC content
#     report inserts with bin assignments for final bin coverage assessment in R
# input:
#     chrom,start0,end1,insert_size_level,frac_gc,n_peak_overlaps on STDIN
#     input sort is not guaranteed, but not required by this script or collate.R
# output:
#     insert size distribution table (n_is_gc) to tmp file
#     indvidual inserts as chr,cb0,isl,peak (ins_bin_isl) to STDOUT

use strict;
use warnings;

# constants
use constant { 
    CHROM             => 0,
    START0            => 1,
    END1              => 2,
    INSERT_SIZE_LEVEL => 3,
    FRAC_GC           => 4,
    N_PEAK_OVERLAPS   => 5,
    #---------------------
    NOT_IN_PEAK    => 0,
    IN_PEAK        => 1
};

# variables
my $TMP_FILE_PREFIX = $ENV{TMP_FILE_PREFIX};
my $FILENAME_PREFIX = $ENV{FILENAME_PREFIX};
my $BIN_SIZE        = $ENV{BIN_SIZE};
my $PRIMARY_GENOME  = $ENV{PRIMARY_GENOME};

# initialize insert size counts
my @insertSizes;

# run the individual inserts
while(<STDIN>){ 
    chomp; 
    my @ins = split("\t"); 

    # parse the genome and chromosome
    my ($chrom, $genome) = split("-", $ins[CHROM]);
    $genome eq $PRIMARY_GENOME or next;

    # increment the appropriate insert size vs. GC table
    my $is = $ins[END1] - $ins[START0];
    my $gcI = int($ins[FRAC_GC] * 100 + 0.5); # convert to integer percent GC
    my $peak = $ins[N_PEAK_OVERLAPS] ? IN_PEAK : NOT_IN_PEAK;
    $insertSizes[$peak][$is][$gcI]++;

    # print individual inserts with bin, insert size level, and counting weight
    print join("\t", 
        $ins[CHROM], # chr + bs0 defins a bin
        int($ins[START0] / $BIN_SIZE) * $BIN_SIZE, # bin start0
        $ins[INSERT_SIZE_LEVEL],
        $peak
    ), "\n";
}

# print the insert size table to temp file
my $isFile = "$TMP_FILE_PREFIX.$FILENAME_PREFIX.n_is_gc.txt";
open my $outH, ">", $isFile or die "recollate_sample_inserts error: cannot open $isFile for writing\n";
foreach my $peak(NOT_IN_PEAK..IN_PEAK){
    foreach my $insertSize($ENV{MIN_INSERT_SIZE}..$ENV{MAX_INSERT_SIZE}){
        if($insertSizes[$peak][$insertSize]){
            print $outH join("\t", $peak, map { $insertSizes[$peak][$insertSize][$_] || 0 } 0..100), "\n";
        } else {
            print $outH join("\t", $peak, (0) x 101), "\n";
        }
    }
}
close $outH;
