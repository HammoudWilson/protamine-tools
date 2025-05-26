# action:
#     for both primary and spike-in genomes:
#       calculate insert size distributions per genome, stratified by percent GC
#       these distributions are used to:
#           assess relative insert sizes of primary and spike-in genomes
#           establish spike-in normalization factors
#     for just the primary genome, where we seek to compare normalized bins to each other:
#       score the actual usage of all possible Tn5 cleavage site kmers to establish sample-level Tn5 site preferences
# input:
#     sorted chrom, start0_padded, end1_padded, insert_size_level, ref_seq on STDIN
# output:
#     insert size distribution table (n_is_gc) to tmp file
#     Tn5 cleavage site kmer counts (n_tn5_smp) to tmp file
#     indvidual inserts as chr,cb0,isl (ins_bin_isl), to STDOUT

use strict;
use warnings;

# constants
use constant { 
    CHROM             => 0,
    START0            => 1,
    END1              => 2,
    INSERT_SIZE_LEVEL => 3,
    REF_SEQ           => 4,
};

# # variables
my $TMP_FILE_PREFIX = $ENV{TMP_FILE_PREFIX};
my $TN5_FLANK_PADDING = $ENV{TN5_FLANK_PADDING};
my @TN5_PREFERENCE_POSITIONS = split(/\s+/, $ENV{TN5_PREFERENCE_POSITIONS});
my $N_TN5_PREFERENCE_POSITIONS = $ENV{N_TN5_PREFERENCE_POSITIONS};
my $TN5_PREFERENCE_SPAN = $ENV{TN5_PREFERENCE_SPAN};
my $nTn5Kmers = 4 ** $N_TN5_PREFERENCE_POSITIONS;
my $BIN_SIZE = $ENV{BIN_SIZE};

# initialize insert size counts and size level file handles
# perform mappability intersection and parsing in these output streams
my %insertSizes;
$insertSizes{$ENV{PRIMARY_GENOME}}  = [];
$insertSizes{$ENV{SPIKE_IN_GENOME}} = [];

# load the possible Tn5 cleavage site kmers
my (%tn5KmerIs, @tn5Kmers);
my @tn5Kmer_counts = (0) x $nTn5Kmers;
open my $inH, "<", $ENV{TN5_KMERS_FILE} or die "collate_count_inserts error: cannot open $ENV{TN5_KMERS_FILE} for reading\n";
my $i = 0;
while(<$inH>){
    chomp;
    push @tn5Kmers, $_;
    $tn5KmerIs{$_} = $i;
    $i++;
}
close $inH;

# run the alignment data
while(<STDIN>){ 
    chomp; 
    my @ins = split("\t"); 

    # unpad insert spans back to the original start0 and end1
    $ins[START0] = $ins[START0] + $TN5_FLANK_PADDING; 
    $ins[END1]   = $ins[END1]   - $TN5_FLANK_PADDING;

    # parse the genome and chromosome
    my ($chrom, $genome) = split("-", $ins[CHROM]);

    # calculate frac GC and increment the insert size vs. GC table
    # do this for both primary and spike-in genomes
    my $gcCount   = $ins[REF_SEQ] =~ tr/GCgc//;
    my $acgtCount = $ins[REF_SEQ] =~ tr/ACGTacgt//;
    $insertSizes{$genome}[$ins[END1] - $ins[START0]][int($gcCount / $acgtCount * 100 + 0.5)]++; 

    # remaining steps apply to primary genome only
    $genome eq $ENV{PRIMARY_GENOME} or next;

    # count cleavage site kmers
    countTn5Kmer(uc substr($ins[REF_SEQ], 0, $TN5_PREFERENCE_SPAN), 0);
    countTn5Kmer(uc substr($ins[REF_SEQ],   -$TN5_PREFERENCE_SPAN), 1);

    # print individual inserts with bin and insert size level metadata
    print join("\t", 
        $ins[CHROM],
        int($ins[START0] / $BIN_SIZE) * $BIN_SIZE, # bin start0
        $ins[INSERT_SIZE_LEVEL]
    ), "\n";
}
sub countTn5Kmer {
    my ($seq, $rc) = @_;
    $rc and rc(\$seq);
    # extract just TN5_PREFERENCE_POSITIONS from $seq
    # thus, the input 13-mer is reduced a non-adjacent 9-mer
    my $kmer = '';
    foreach my $pos (@TN5_PREFERENCE_POSITIONS) {
        $kmer .= substr($seq, $pos, 1);
    }
    exists $tn5KmerIs{$kmer} or return;
    $tn5Kmer_counts[$tn5KmerIs{$kmer}]++;
}
sub rc { 
    my ($seq) = @_;
    $$seq = reverse $$seq;
    $$seq =~ tr/ATGC/TACG/;
}

# print the insert size table to temp file
my $isFile = "$TMP_FILE_PREFIX.n_is_gc.txt";
open my $outH, ">", $isFile or die "collate_count_inserts error: cannot open $isFile for writing\n";
foreach my $genome($ENV{PRIMARY_GENOME}, $ENV{SPIKE_IN_GENOME}){ 
    foreach my $insertSize($ENV{MIN_INSERT_SIZE}..$ENV{MAX_INSERT_SIZE}){
        if($insertSizes{$genome} and $insertSizes{$genome}[$insertSize]){
            print $outH join("\t", $genome, map { $insertSizes{$genome}[$insertSize][$_] || 0 } 0..100), "\n";
        } else {
            print $outH join("\t", $genome, (0) x 101), "\n";
        }
    }
}
close $outH;

# print the Tn5 cleavage site kmers to temp BED file
my $tn5File = "$TMP_FILE_PREFIX.n_tn5_smp.txt";
open $outH, ">", $tn5File or die "collate_count_inserts error: cannot open $tn5File for writing\n";
foreach my $kmer(@tn5Kmers){
    print $outH join("\t", $kmer,  $tn5Kmer_counts[$tn5KmerIs{$kmer}]), "\n";
}
close $outH;

# open $outH, ">", $tn5File or die "collate_count_inserts error: cannot open $tn5File for writing\n";
# foreach my $chrom(@chroms) {
#     foreach my $binI0(0..$maxBinI0{$chrom}) {
#         print $outH join("\t", 
#             $chrom, 
#             $binI0 * $BIN_SIZE, 
#             ($binI0 + 1) * $BIN_SIZE, 
#             $binCounts{$chrom}[$binI0]
#         ), "\n";
#     }
# }
# close $outH;
