# action:
#     step through a masked genome fasta stream
#     extract the Tn5 cleavage site nonamers possible start positions
#     keep and count known nonamers within $TN5_KMERS_FILE
# input:
#     genome FASTA stream on STDIN with all never-mappable bases masked as N
# output:
#     count of Tn5 cleavage site usage on STDOUT

use strict;
use warnings;

# dependencies
require "$ENV{MODULES_DIR}/tn5/tn5_sites.pl";
use vars qw(@insertSizeLevelIs);

# variables
my $TN5_PREFERENCE_SPAN = $ENV{TN5_PREFERENCE_SPAN};
my $INSERT_SIZE_LEVEL = $ENV{INSERT_SIZE_LEVEL};
my $islI = $insertSizeLevelIs[$INSERT_SIZE_LEVEL];
my ($seqBuffer, $pos, $chrom, $genome, $allowedChrom) = ('', 0);

# constants
use constant MAX_BUFFER_POS => 1e6;

# run the alignment data
while(my $line = <STDIN>){ 
    chomp $line; 

    # reset on entering a new chromosome
    if($line =~ /^>/) { 
        $line =~ s/^>//;
        ($chrom, $genome) = split("-", $line);
        $seqBuffer = '';
        $pos = 0;
        $allowedChrom = (
            $chrom !~ m/_/ and
            $chrom ne "chrM" and
            $chrom ne "chrEBV"
        );

    # append a new sequence chunk and process all available candidate Tn5 cleavage sites
    } else {
        $allowedChrom or next;
        $seqBuffer .= $line;
        while(length($seqBuffer) - $pos >= $TN5_PREFERENCE_SPAN) {
            countTn5Site($genome, $islI, uc substr($seqBuffer, $pos, $TN5_PREFERENCE_SPAN));
            $pos++;
            if ($pos > MAX_BUFFER_POS) {
                $seqBuffer = substr($seqBuffer, $pos);
                $pos = 0;
            }
        }
    }
}

# print the aggregated site counts
printTn5SiteCounts();
