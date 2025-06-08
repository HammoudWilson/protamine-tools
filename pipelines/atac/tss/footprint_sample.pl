# action:
#     filter inserts to primary genome only
#     remove genome suffix on chromosome names
#     parse inserts with insert size length and Tn5 kmers into endpoint site weights
# input:
#     chrom,start0,end1,insert_size_level,tn5_site_left,tn5_site_right on STDIN
#     input sort is not guaranteed, but not required by this script (sort occurs after this script)
# outputs:
#     chrom,start0,end1,tn5_weight_left,tn5_weight_right on STDOUT

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
};

# variables
my $PRIMARY_GENOME = $ENV{PRIMARY_GENOME};
my @MAPPABILITY_SIZE_LEVELS = split(/\s+/, $ENV{MAPPABILITY_SIZE_LEVELS});
my ($nInserts, $wInserts) = (0, 0); # counters for inserts and weighted inserts

# initialize insert size level indices
my (@insertSizeLevelIs);
my $islI = 0;
foreach my $insertSizeLevel(@MAPPABILITY_SIZE_LEVELS){
    $insertSizeLevelIs[$insertSizeLevel] = $islI;
    $islI++;
}

# load the possible Tn5 cleavage site kmers
my (%tn5KmerIs);
open my $inH, "<", $ENV{TN5_KMERS_FILE} or die "collate_sample_inserts error: cannot open $ENV{TN5_KMERS_FILE} for reading\n";
my $siteI = 0;
while(<$inH>){
    chomp;
    $tn5KmerIs{$_} = $siteI;
    $siteI++;
}
close $inH;

# load the site frequencies by reference genome when encountered in insert stream
my (@siteFreqObs, @siteFreqExp);
sub loadTn5SiteFreqsFile {
    my ($file, $freqs) = @_;
    open my $fh, "-|", "zcat $file" or die "collate_sample_inserts error: cannot open $file for reading\n";
    while(<$fh>){
        chomp;
        push @$freqs, [split(",")]; # yes, it is comma-separated
    }
    close $fh;
}
my $obsFreqsFile = "$ENV{DATA_FILE_PREFIX}.primary.tn5_site_freqs_obs.txt.gz";
my $expFreqsFile = "$ENV{GENOME_METADATA_PREFIX}.primary.tn5_site_freqs_exp.txt.gz";
loadTn5SiteFreqsFile($obsFreqsFile, \@siteFreqObs);
loadTn5SiteFreqsFile($expFreqsFile, \@siteFreqExp);

# run the individual inserts
while(<STDIN>){ 
    chomp; 
    my @ins = split("\t"); 

    # parse the genome and chromosome
    my ($chrom, $genome) = split("-", $ins[CHROM]);
    $genome eq $PRIMARY_GENOME or next;

    # calculate the insert Tn5 site weights
    #   divide exp / obs to weight less Tn5-preferred sites more heavily when summing position counts
    #   "expected" refers to the Tn5 cleavage site frequency in the genome
    #   "observed" refers to sites found in the sample data)
    #   neither can be zero since tn5LI was a kmer found in the genome, but ratio can be unstable
    #   instability limits not enforced here, this is deferred to handling in the app
    my $tn5LI = $tn5KmerIs{$ins[TN5_SITE_LEFT]};
    my $tn5RI = $tn5KmerIs{$ins[TN5_SITE_RIGHT]};
    (!defined $tn5LI or !defined $tn5RI) and next; # only count inserts with known Tn5 sites, i.e., no N bases allowed at ends
    my $islI = $insertSizeLevelIs[$ins[INSERT_SIZE_LEVEL]];
    my $wStart = int($siteFreqExp[$tn5LI][$islI] / $siteFreqObs[$tn5LI][$islI] * 1000 + 0.5) / 1000;
    my $wEnd   = int($siteFreqExp[$tn5RI][$islI] / $siteFreqObs[$tn5RI][$islI] * 1000 + 0.5) / 1000;

    # print individual inserts with bin and endpoint weights
    print join("\t", 
        $chrom, # simplify output chrom names
        $ins[START0],
        $ins[END1],
        $wStart,
        $wEnd
    ), "\n";
    $nInserts++;
    $wInserts += ($wStart + $wEnd) / 2;
}

# print tallies for normalizing visualization plot
open my $outH, ">", $ENV{TALLY_FILE} or die "footprint_samples error: cannot open $ENV{TALLY_FILE} for writing\n";
print $outH join("\t", $ENV{BGZ_FILE_NAME}, $nInserts, $wInserts), "\n";
close $outH;
