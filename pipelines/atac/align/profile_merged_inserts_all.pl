
# this script is not executed within a pipeline
# it is used in BAM_OUTPUT_DIR as
#     samtools cat *.bam | samtools view - | perl profile_merged_inserts_all.pl

# action:
#   extract an insert-size-dependent mappability profile of merged reads
#   used to assess mappability and GC content as a function of discrete insert size levels
#   works on merged reads only, which encompass mono- and sub-nucleosomal reads
# input:
#     SAM stream of Bowtie2 alignments on STDIN
# outputs:
#     insert size distribution table written to STDOUT

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
    _SKIP_ALIGNMENT => 3841, # _PAIRED + _NOT_PRIMARY + _FAILED_QC + _DUPLICATE + _SUPPLEMENTARY
};

# variables
my $MIN_MAPQ = 30;
my @MAPPABILITY_KMER_LENGTHS = split(/\s+/, "35 40 45 50 55 60 70 80 90 100 120 140 160 180 200 220 240 260 280 300");
my %counts = map {
    $_ => {
        n_inserts                  => 0,
        n_inserts_unaligned        => 0,
        n_inserts_aligned_failed   => 0,
        n_inserts_aligned_passed   => 0,
        gc_bases_unaligned         => 0,
        gc_bases_aligned_failed    => 0,
        gc_bases_aligned_passed    => 0,
        acgt_bases_unaligned       => 0,
        acgt_bases_aligned_failed  => 0,
        acgt_bases_aligned_passed  => 0,
    }
 } @MAPPABILITY_KMER_LENGTHS;

# run the alignment data
while(my $line = <STDIN>){ 
    my @aln = split("\t", $line, SPLIT_TO_SEQ);
    ($aln[FLAG] & _SKIP_ALIGNMENT) and next; 

    my $insertSizeLevel = getInsertSizeLevel(length($aln[SEQ])) or next; # ignore teeny-tiny inserts
    $counts{$insertSizeLevel}{n_inserts}++;

    my $gcCount   = $aln[SEQ] =~ tr/GCgc//;
    my $acgtCount = $aln[SEQ] =~ tr/ACGTacgt//;

    if($aln[FLAG] & _UNMAPPED){
        $counts{$insertSizeLevel}{n_inserts_unaligned}++;
        $counts{$insertSizeLevel}{gc_bases_unaligned}   += $gcCount;
        $counts{$insertSizeLevel}{acgt_bases_unaligned} += $acgtCount;
    } elsif($aln[MAPQ] < $MIN_MAPQ){
        $counts{$insertSizeLevel}{n_inserts_aligned_failed}++;
        $counts{$insertSizeLevel}{gc_bases_aligned_failed}   += $gcCount;
        $counts{$insertSizeLevel}{acgt_bases_aligned_failed} += $acgtCount;
    } else {
        $counts{$insertSizeLevel}{n_inserts_aligned_passed}++;
        $counts{$insertSizeLevel}{gc_bases_aligned_passed}   += $gcCount;
        $counts{$insertSizeLevel}{acgt_bases_aligned_passed} += $acgtCount;
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
foreach my $insert_size_level (@MAPPABILITY_KMER_LENGTHS) {
    my $profile = $counts{$insert_size_level};

    my $frac_unaligned      = $$profile{n_inserts} ? $$profile{n_inserts_unaligned}      / $$profile{n_inserts} : 0;
    my $frac_aligned_failed = $$profile{n_inserts} ? $$profile{n_inserts_aligned_failed} / $$profile{n_inserts} : 0;
    my $frac_aligned_passed = $$profile{n_inserts} ? $$profile{n_inserts_aligned_passed} / $$profile{n_inserts} : 0;

    my $frac_gc_unaligned      = $$profile{acgt_bases_unaligned}      ? $$profile{gc_bases_unaligned}      / $$profile{acgt_bases_unaligned}      : 0;
    my $frac_gc_aligned_failed = $$profile{acgt_bases_aligned_failed} ? $$profile{gc_bases_aligned_failed} / $$profile{acgt_bases_aligned_failed} : 0;
    my $frac_gc_aligned_passed = $$profile{acgt_bases_aligned_passed} ? $$profile{gc_bases_aligned_passed} / $$profile{acgt_bases_aligned_passed} : 0;

    print join("\t", 
        $insert_size_level,

        $$profile{n_inserts},
        $$profile{n_inserts_unaligned},
        $$profile{n_inserts_aligned_failed},
        $$profile{n_inserts_aligned_passed},

        sprintf("%.3f", $frac_unaligned),
        sprintf("%.3f", $frac_aligned_failed),
        sprintf("%.3f", $frac_aligned_passed),

        sprintf("%.3f", $frac_gc_unaligned),
        sprintf("%.3f", $frac_gc_aligned_failed),
        sprintf("%.3f", $frac_gc_aligned_passed),
    ), "\n";
}
