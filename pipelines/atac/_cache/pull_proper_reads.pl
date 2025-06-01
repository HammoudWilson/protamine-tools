# action:
#     restrict to proper read pairs with no SVs (whether merged or not)
# input:
#     coordinate-sorted SAM stream with header on STDIN
# outputs:
#     filtered coordinate-sorted SAM stream with header on STDOUT
#       one read per line with positive TLEN

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

# run the alignment data
while(my $line = <STDIN>){ 

    # pass header lines as-is
    if($line =~ /^@/){
        print $line;
        next;
    }

    # process alignments by QNAME
    # remember, input stream is coordinate sorted, cannot group by QNAME
    my @aln = split("\t", $line, SPLIT_TO_TLEN); 
    ($aln[FLAG] & _SUPPLEMENTARY) and next; # redundant since Bowtie2 used in end-to-end mode
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

    # print just one filtered and/or updated alignment per read to STDOUT
    print join("\t", @aln); 
}
