# action:
#   extract an insert-size-dependent mappability profile of merged reads
#   used to assess mappability and GC content as a function of discrete insert size levels
#   applied to merged reads only, which encompass mono- and sub-nucleosomal reads
# input:
#     SAM stream of Bowtie2 alignments of merged reads on STDIN
# outputs:
#     SAM stream repeated verbatim on STDOUT
#     insert size distribution table written to sample-level $MAPPABILITY_PROFILE file
#       $FILENAME_PREFIX,insert_size_level,total_inserts,mappable_inserts,frac_mappable,frac_GC_mappable,frac_GC_unmappable
#       where:
#           $FILENAME_PREFIX allows concatenation of multiple samples
#           insert_size_levels taken from $MAPPABILITY_KMER_LENGTHS
#           mappable_inserts have MAPQ >= $MIN_MAPQ (not FLAG & _UNMAPPED, see note below)
#           frac_mappable is mappable_inserts divided by total_inserts
#           frac_GC_xxx is the total GC base count divided by the total base count per insert_size_level (same as mean(frac_GC))

use strict;
use warnings;

# constants
use constant { 
    QNAME   => 0,
    FLAG    => 1,
    RNAME   => 2,
    POS     => 3,
    MAPQ    => 4,
    CIGAR   => 5,
    RNEXT   => 6,
    PNEXT   => 7,
    TLEN    => 8,
    SEQ     => 9,
    QUAL    => 10,
    TAGS    => 11,
    #-----------------
    SPLIT_TO_SEQ => 11,
    #-----------------
    _PAIRED         => 0x1,
    _PROPER_PAIR    => 0x2,
    _UNMAPPED       => 0x4,
    _MATE_UNMAPPED  => 0x8,
    _REVERSE        => 0x10,
    _MATE_REVERSE   => 0x20,
    _FIRST_IN_PAIR  => 0x40,
    _SECOND_IN_PAIR => 0x80,
    _NOT_PRIMARY    => 0x100,
    _FAILED_QC      => 0x200,
    _DUPLICATE      => 0x400,
    _SUPPLEMENTARY  => 0x800,
    #-----------------
    _SKIP_ALIGNMENT => 1796, # _UNMAPPED + _NOT_PRIMARY + _FAILED_QC + _DUPLICATE
};

# variables
my ($prevQname, @alns);
my $MIN_MAPQ = $ENV{MIN_MAPQ};
my @MAPPABILITY_KMER_LENGTHS = split(/\s+/, $ENV{MAPPABILITY_KMER_LENGTHS});
my %counts = map {
    $_ => {
        total_inserts           => 0,
        mappable_inserts        => 0,
        gc_bases_mappable       => 0,
        gc_bases_unmappable     => 0,
        acgt_bases_mappable     => 0,
        acgt_bases_unmappable   => 0,
    }
 } @MAPPABILITY_KMER_LENGTHS;

# run the alignment data
while(my $line = <STDIN>){ 
    print $line; # repeat all lines to STDOUT
    next if $line =~ /^@/; # skip header lines
    my @aln = split("\t", $line, SPLIT_TO_SEQ);

    # skip _UNMAPPED reads
    #   they are not "unmappable" in mappability context (those are alignments with low MAPQ that cannot be uniquely placed)
    #   instead, _UNMAPPED are "unusable", e.g., low quality sequences that could not be recognized as anything
    # skip if any failure flag is set
    # !_PAIRED enforced by placement of this script in a stream
    #   _PROPER_PAIR, _MATE_UNMAPPED, _MATE_REVERSE, _FIRST_IN_PAIR, _SECOND_IN_PAIR therefore not relevant
    # _SUPPLEMENTARY handled by alignment counting below
    $aln[FLAG] & _SKIP_ALIGNMENT and next; 

    # process read alignment groups
    if($prevQname and $prevQname ne $aln[QNAME]){
        processQname();
        @alns = ();
    }
    push @alns, \@aln;
    $prevQname = $aln[QNAME];
}
# process the last read alignment group
processQname();

sub processQname {
    @alns == 1 or return; # read must not have SVs to be profiled
    my $insertSizeLevel = getInsertSizeLevel(length($alns[0][SEQ])) or return; # ignore teeny-tiny inserts
    $counts{$insertSizeLevel}{total_inserts}++;
    my $gcCount   = $alns[0][SEQ] =~ tr/GCgc//;
    my $acgtCount = $alns[0][SEQ] =~ tr/ACGTacgt//;
    if($alns[0][MAPQ] >= $MIN_MAPQ){
        $counts{$insertSizeLevel}{mappable_inserts}++;
        $counts{$insertSizeLevel}{gc_bases_mappable}   += $gcCount;
        $counts{$insertSizeLevel}{acgt_bases_mappable} += $acgtCount;
    } else { # for this script, unmappable means low MAPQ, not alignment failure
        $counts{$insertSizeLevel}{gc_bases_unmappable}   += $gcCount;
        $counts{$insertSizeLevel}{acgt_bases_unmappable} += $acgtCount;
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

# print the mappability profile
open my $outH, ">", $ENV{MAPPABILITY_PROFILE} or die "profile_merged_inserts error: cannot open $ENV{MAPPABILITY_PROFILE} for writing\n";
foreach my $insert_size_level (@MAPPABILITY_KMER_LENGTHS) {
    my $profile = $counts{$insert_size_level};
    my $frac_mappable = $$profile{total_inserts} ? $$profile{mappable_inserts} / $$profile{total_inserts} : 0;
    my $frac_GC_mappable   = $$profile{acgt_bases_mappable}   ? $$profile{gc_bases_mappable}   / $$profile{acgt_bases_mappable}   : 0;
    my $frac_GC_unmappable = $$profile{acgt_bases_unmappable} ? $$profile{gc_bases_unmappable} / $$profile{acgt_bases_unmappable} : 0;
    print $outH join("\t", 
        $ENV{FILENAME_PREFIX},
        $insert_size_level,
        $$profile{total_inserts},
        $$profile{mappable_inserts},
        sprintf("%.3f", $frac_mappable),
        sprintf("%.3f", $frac_GC_mappable),
        sprintf("%.3f", $frac_GC_unmappable)
    ), "\n";
}
close $outH;
