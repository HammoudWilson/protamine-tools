# action:
#     account for Tn5 cleavage site preference by calculating counting weights
#     calculate insert size distributions per genome and site weighting, stratified by insert percent GC
#     these distributions are used to:
#       assess relative insert sizes of primary and spike-in genomes
#       establish spike-in normalization factors
#       explore the correlation of insert size and GC content, before and after Tn5 cleavage site adjustment
#     report inserts with bin assignments and site weights for final bin coverage assessment in R
# input:
#     chrom,start0,end1,insert_size_level,tn5_site_left,tn5_site_right,frac_gc on STDIN
#     input sort is not guaranteed, but not required by this script or collate.R
# output:
#     insert size distribution table (n_is_gc) to tmp file
#     indvidual inserts as chr,cb0,isl,wgt (ins_bin_isl) to STDOUT

use strict;
use warnings;

# constants
use constant { 
    CHROM             => 0,
    START0            => 1,
    END1              => 2,
    INSERT_SIZE_LEVEL => 3,
    TN5_SITE_LEFT     => 4,
    TN5_SITE_RIGHT    => 5,
    FRAC_GC           => 6,
    #----------------------
    OBSERVED => 0, # weighting types used when tabulating insert size by GC
    WEIGHTED => 1
};

# variables
my $TMP_FILE_PREFIX = $ENV{TMP_FILE_PREFIX};
my $FILENAME_PREFIX = $ENV{FILENAME_PREFIX};
my $BIN_SIZE = $ENV{BIN_SIZE};
my @MAPPABILITY_SIZE_LEVELS = split(/\s+/, $ENV{MAPPABILITY_SIZE_LEVELS});

# initialize insert size counts
my %insertSizes;
$insertSizes{$ENV{PRIMARY_GENOME} }[OBSERVED] = [];
$insertSizes{$ENV{PRIMARY_GENOME} }[WEIGHTED] = [];
$insertSizes{$ENV{SPIKE_IN_GENOME}}[OBSERVED] = [];
$insertSizes{$ENV{SPIKE_IN_GENOME}}[WEIGHTED] = [];

# initialize insert size level indices
my (@insertSizeLevelIs);
my $islI = 0;
foreach my $insertSizeLevel(@MAPPABILITY_SIZE_LEVELS){
    $insertSizeLevelIs[$insertSizeLevel] = $islI;
    $islI++;
}

# load the possible Tn5 cleavage site kmers
my (@tn5Kmers, %tn5KmerIs);
open my $inH, "<", $ENV{TN5_KMERS_FILE} or die "collate_sample_inserts error: cannot open $ENV{TN5_KMERS_FILE} for reading\n";
my $siteI = 0;
while(<$inH>){
    chomp;
    push @tn5Kmers, $_;
    $tn5KmerIs{$_} = $siteI;
    $siteI++;
}
close $inH;

# load the site frequencies by reference genome when encountered in insert stream
my ($workingGenome, @siteFreqObs, @siteFreqExp) = ("");
sub loadTn5SiteFreqsFile {
    my ($file, $freqs) = @_;
    open my $fh, "-|", "zcat $file" or die "collate_sample_inserts error: cannot open $file for reading\n";
    while(<$fh>){
        chomp;
        push @$freqs, [split(",")]; # yes, it is comma-separated
    }
    close $fh;
}
sub loadTn5SiteFreqs {
    my ($genome) = @_;
    $workingGenome = $genome;
    my $refType = $genome eq $ENV{PRIMARY_GENOME} ? "primary" : "spike_in";
    my $obsFreqsFile = "$ENV{DATA_FILE_PREFIX}.$refType.tn5_site_freqs_obs.txt.gz";
    my $expFreqsFile = "$ENV{GENOME_METADATA_PREFIX}.$refType.tn5_site_freqs_exp.txt.gz";
    @siteFreqObs = @siteFreqExp = ();
    loadTn5SiteFreqsFile($obsFreqsFile, \@siteFreqObs);
    loadTn5SiteFreqsFile($expFreqsFile, \@siteFreqExp);
}

# run the individual inserts
while(<STDIN>){ 
    chomp; 
    my @ins = split("\t"); 

    # parse the genome and chromosome
    my ($chrom, $genome) = split("-", $ins[CHROM]);
    $genome eq $workingGenome or loadTn5SiteFreqs($genome);

    # calculate the summable insert Tn5 site weights
    #   multiply site frequencies at each independent insert end
    #   divided exp / obs to weight less preferred sites more heavily in bin counts
    my $tn5LI = $tn5KmerIs{$ins[TN5_SITE_LEFT]};
    my $tn5RI = $tn5KmerIs{$ins[TN5_SITE_RIGHT]};
    (!defined $tn5LI or !defined $tn5RI) and next; # only count inserts with known Tn5 sites, i.e., no N bases allowed at ends
    my $islI = $insertSizeLevelIs[$ins[INSERT_SIZE_LEVEL]];
    my $weight = ($siteFreqExp[$tn5LI][$islI] * $siteFreqExp[$tn5RI][$islI]) / 
                 ($siteFreqObs[$tn5LI][$islI] * $siteFreqObs[$tn5RI][$islI]);  

    # increment the appropriate insert size vs. GC table
    my $is = $ins[END1] - $ins[START0];
    my $gcI = int($ins[FRAC_GC] * 100 + 0.5); # convert to integer percent GC
    $insertSizes{$genome}[OBSERVED][$is][$gcI]++;
    $insertSizes{$genome}[WEIGHTED][$is][$gcI] += $weight;

    # print individual inserts with bin, insert size level, and counting weight
    print join("\t", 
        $ins[CHROM], # chr + bs0 defins a bin
        int($ins[START0] / $BIN_SIZE) * $BIN_SIZE, # bin start0
        $ins[INSERT_SIZE_LEVEL],
        $weight # Tn5 cleavage site adjusted counting weight
    ), "\n";
}

# print the insert size table to temp file
my $isFile = "$TMP_FILE_PREFIX.$FILENAME_PREFIX.n_is_gc.txt";
open my $outH, ">", $isFile or die "collate_sample_inserts error: cannot open $isFile for writing\n";
foreach my $genome($ENV{PRIMARY_GENOME}, $ENV{SPIKE_IN_GENOME}){
    foreach my $weighting(OBSERVED, WEIGHTED){
        my $wgtName = $weighting == OBSERVED ? "observed" : "weighted";
        my $is = $insertSizes{$genome}[$weighting];
        foreach my $insertSize($ENV{MIN_INSERT_SIZE}..$ENV{MAX_INSERT_SIZE}){
            if($$is[$insertSize]){
                print $outH join("\t", $genome, $wgtName, map { $$is[$insertSize][$_] || 0 } 0..100), "\n";
            } else {
                print $outH join("\t", $genome, $wgtName, (0) x 101), "\n";
            }
        }
    }
}
close $outH;
