# action:
#     score the actual usage of all possible Tn5 cleavage site nonamers to establish sample-level Tn5 site preferences
#     do this for both the primary and spike-in genomes to assess a spermatid-specific impact on Tn5 sites
#     do this stratified by insert size level, to account for different mappability spans at different insert size levels
#     write a file with insert spans per sample for use by collate
# input:
#     sorted chrom, start0_padded, end1_padded, insert_size_level, ref_seq on STDIN
# output:
#     insert spans file in BED format, compressed with gzip, with Tn5 cleavage site nonamers
#     Tn5 cleavage site nonamer counts by insert size level to STDOUT
#       one row per nonamer sorted by genome and nonamer, columns sorted by insert size level

use strict;
use warnings;

# constants
use constant { 
    CHROM             => 0,
    START0_PADDED     => 1, # coordinates and REF_SEQ are padded with TN5_FLANK_PADDING bases
    END1_PADDED       => 2,
    INSERT_SIZE_LEVEL => 3,
    REF_SEQ_PADDED    => 4,
};

# variables
my $TN5_PREFERENCE_SPAN = $ENV{TN5_PREFERENCE_SPAN};
my $TN5_FLANK_PADDING = $ENV{TN5_FLANK_PADDING};

# insert spans output file
# despite the upstream sort merge, the final ins_dedup_mpp.bed.gz insert_spans file is not guaranteed to be sorted (and is not)
#   unclear why, possibly bedtools getfasta does this?
#   at present, no downstream code requires sorted insert spans but caution should be taken
my $outFile = "$ENV{INSERT_SPANS_DIR}/$ENV{FILENAME_PREFIX}.ins_dedup_mpp.bed.gz";
open my $insH, "|-", "gzip -c > $outFile" or 
    die "count_cleavage_sites error: cannot open $outFile for writing\n";

# dependencies
require "$ENV{MODULES_DIR}/tn5/tn5_sites.pl";
use vars qw(@insertSizeLevelIs);

# run the alignment data
while(<STDIN>){ 
    chomp; 
    my @ins = split("\t"); 
    my ($chrom, $genome) = split("-", $ins[CHROM]);
    my $islI = $insertSizeLevelIs[$ins[INSERT_SIZE_LEVEL]];
    my $site_left  = countTn5Site($genome, $islI, uc substr($ins[REF_SEQ_PADDED], 0, $TN5_PREFERENCE_SPAN)) or next;
    my $site_right = countTn5Site($genome, $islI, uc substr($ins[REF_SEQ_PADDED],   -$TN5_PREFERENCE_SPAN)) or next;
    print $insH join("\t",
        $ins[CHROM], 
        $ins[START0_PADDED] + $TN5_FLANK_PADDING, # unpad coordinates back to original insert span on reference
        $ins[END1_PADDED]   - $TN5_FLANK_PADDING, 
        $ins[INSERT_SIZE_LEVEL], # insert size level and cleavage site used by atac/collate to weight inserts when counting
        $site_left,
        $site_right
    ), "\n";
}
close $insH;

# print the aggregated site counts
printTn5SiteCounts();
