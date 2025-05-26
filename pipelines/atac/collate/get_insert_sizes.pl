# action:
#     unpack the distinct insert end1s for every insert start0
#     calculate insert size distributions per genome, stratified by percent GC
# input:
#     non-yet-deduplicated BED4 on STDIN with fraction GC in 4th column
# outputs:
#     insert size distribution table on STDOUT

use strict;
use warnings;

# constants
use constant { 
    CHROM   => 0,
    START0  => 1,
    END1    => 2,
    FRAC_GC => 3,
};

# variables
my %counts;
$counts{$ENV{PRIMARY_GENOME}}  = [];
$counts{$ENV{SPIKE_IN_GENOME}} = [];

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
# process the last START0
processStart0();

# dupllicate and count insert sizes in a set of alignments starting at the same START0
sub processStart0 {
    my %end1s;
    foreach my $f(@start0Group){
        $end1s{$$f[END1]} and next; # deduplication occurs here
        $end1s{$$f[END1]}++;
        my ($chrom, $genome) = split("-", $$f[CHROM]);
        $counts{$genome}[$$f[END1] - $$f[START0]][int($$f[FRAC_GC] * 100 + 0.5)]++; 
    }
}

# print the counts
foreach my $genome($ENV{PRIMARY_GENOME}, $ENV{SPIKE_IN_GENOME}){ 
    foreach my $insertSize($ENV{MIN_INSERT_SIZE}..$ENV{MAX_INSERT_SIZE}){
        if($counts{$genome} and $counts{$genome}[$insertSize]){
            print join("\t", $genome, map { $counts{$genome}[$insertSize][$_] || 0 } 0..100), "\n";
        } else {
            print join("\t", $genome, (0) x 101), "\n";
        }
    }
}
