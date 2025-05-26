# action:
#     unpack and de-duplicate the distinct insert end1s for every insert start0
#     perform size-level-dependent intersection to genmap mappability
# input:
#     non-yet-deduplicated BED3 on STDIN
# outputs:
#     deduplicated inserts, intersected against mappability, to tmp files stratified by insert size level

use strict;
use warnings;

# constants
use constant { 
    CHROM   => 0,
    START0  => 1,
    END1    => 2,
};

# variables
my @MAPPABILITY_KMER_LENGTHS = split(/\s+/, $ENV{MAPPABILITY_KMER_LENGTHS});
my $nPassed = 0;

# initialize insert size level file handles
# perform mappability intersection and parsing in these output streams
my (%insertSizeLevelHs);
foreach my $insertSizeLevel(@MAPPABILITY_KMER_LENGTHS){
    my $islFile = "$ENV{TMP_FILE_PREFIX}.ins_dedup_mpp_$insertSizeLevel.txt";
    my ($mppFile) = glob("$ENV{MAPPABILITY_FILE_PREFIX}.k_$insertSizeLevel.e_*.bed.gz");
    my $intersect = "bedtools intersect -u -sorted -g $ENV{GENOME_FASTA_SHM}.fai -a - -b $mppFile | awk 'BEGIN{OFS=\"\\t\"}{print \$1, \$2, \$4, \$5}' > $islFile";
    open my $outH, "|-", $intersect or die "collate_dedup_mappable error: cannot open $islFile for streaming\n";
    $insertSizeLevelHs{$insertSizeLevel} = $outH;
}

# run the alignment data
my ($prevKey, @start0Group);
while(<STDIN>){ 
    chomp; 
    my @f = split("\t"); 
    my $key = join("\t", @f[CHROM, START0]); # input bam is pre-sorted by chrom and start0
    if($prevKey and $prevKey ne $key){
        processStart0();
        @start0Group = ();
    }
    push @start0Group, \@f;
    $prevKey = $key;
}
processStart0();

# process alignments starting at the same START0
sub processStart0 {
    my %end1s;
    foreach my $f(@start0Group){

        # deduplicate inserts to a count of one insert per unique genome span
        $end1s{$$f[END1]} and next; 
        $end1s{$$f[END1]}++;
        $nPassed++;

        # convert insert size to insert size level
        my $insertSizeLevel = getInsertSizeLevel($$f[END1] - $$f[START0]) or next;

        # print the bin counts to the proper insert size level file
        my $outH = $insertSizeLevelHs{$insertSizeLevel};
        print $outH join("\t", 
            $$f[CHROM],
            $$f[START0],
            $$f[START0] + 1, # thus each left-justified insert start can now be interected with genmap mappability ...
            $$f[END1],       # ... while retaining each original insert end (could be multiple for each start)
            $insertSizeLevel
        ), "\n";
    }
}
sub getInsertSizeLevel {
    my ($insertSize) = @_;
    my $level = 0;
    foreach my $kmerLength (@MAPPABILITY_KMER_LENGTHS) {
        $kmerLength > $insertSize and last;
        $level = $kmerLength;
    }
    return $level;
}

# close all insert ouput handles
foreach my $insertSizeLevel(@MAPPABILITY_KMER_LENGTHS){
    my $outH = $insertSizeLevelHs{$insertSizeLevel};
    close $outH;
}

# print the dedup count to stage count file
open my $stageCountH, ">>", $ENV{STAGE_COUNT_FILE} or die "collate_dedup_mappable error: cannot open $ENV{STAGE_COUNT_FILE} for appending\n";
print $stageCountH join("\t", $ENV{SAMPLE_LABEL}, "stage_3_dedup", $nPassed), "\n";
close $stageCountH;
