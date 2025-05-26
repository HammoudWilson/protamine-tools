# action:
#     remove SAM alignments on non-nuclear chroms
#     restrict to proper read pairs with no SVs (whether merged or not)
#     identify the genome span of the (merged) read pair
#     reject reads with genomic spans outside the expected insert size limits
# input:
#     coordinate-sorted SAM stream on STDIN
# outputs:
#     filtered, coordinate-sorted BED3 stream on STDOUT
#       one line per (merged) read pair span, with exact reference coordinates

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
    #------------------
    SPLIT_TO_TLEN => 10,
    #------------------
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
};

# variables
my $MIN_INSERT_SIZE = $ENV{MIN_INSERT_SIZE};
my $MAX_INSERT_SIZE = $ENV{MAX_INSERT_SIZE}; 
my %allowed;
foreach my $chrom(1..22){
    $allowed{"chr$chrom-$ENV{PRIMARY_GENOME}"}  = 1; # composite chrom names
    $allowed{"chr$chrom-$ENV{SPIKE_IN_GENOME}"} = 1;
}
foreach my $chrom(qw(X Y 2L 2R 3L 3R)){
    $allowed{"chr$chrom-$ENV{PRIMARY_GENOME}"}  = 1;
    $allowed{"chr$chrom-$ENV{SPIKE_IN_GENOME}"} = 1;
}

# run the alignment data
# remember, input stream is coordinate sorted, cannot group by QNAME
while(my $line = <STDIN>){ 
    my @aln = split("\t", $line, SPLIT_TO_TLEN); 
    $allowed{$aln[RNAME]} or next;
    ($aln[FLAG] & _SUPPLEMENTARY) and next; # redundant, since Bowtie2 used in end-to-end mode
    if($aln[FLAG] & _PAIRED){
        ($aln[FLAG] & _PROPER_PAIR) or next;
        ($aln[FLAG] & _REVERSE) and next; # thus, all output TLEN are positive
    } else {
        $aln[TLEN] = do { # calculate TLEN for merged reads
            my $nRefBases = 0;
            while ($aln[CIGAR] =~ /(\d+)([MD])/g) {
                $nRefBases += $1;
            }
            $nRefBases;
        };
    }

    # enforce insert size limits
    # redundant, these should always pass given fastp and bowtie2 handling upstream
    $aln[TLEN] < $MIN_INSERT_SIZE and next;
    $aln[TLEN] > $MAX_INSERT_SIZE and next;

    # print the read pair genome span as BED3 to STDOUT
    print join("\t", 
        $aln[RNAME], 
        $aln[POS] - 1, 
        $aln[POS] + $aln[TLEN] - 1
    ), "\n"; 
}
