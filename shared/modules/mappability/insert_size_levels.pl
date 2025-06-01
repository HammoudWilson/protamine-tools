# support for assessing mappability insert size levels

use strict;
use warnings;

# variables 
my @MAPPABILITY_SIZE_LEVELS = split(/\s+/, $ENV{MAPPABILITY_SIZE_LEVELS});

# return the insert size level as the floor of the insert size
# i.e., the largest size level less than or equal to the insert size
sub getInsertSizeLevel {
    my ($insertSize) = @_;
    my $level = 0;
    foreach my $kmerLength (@MAPPABILITY_SIZE_LEVELS) {
        $kmerLength > $insertSize and last;
        $level = $kmerLength;
    }
    return $level;
}

1;
