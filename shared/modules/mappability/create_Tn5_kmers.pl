# action:
#     create a file with all possible Tn5 kmers
#     these are 9-mers representing non-adjacent positions relative to Tn5 cleavage positions
#       see set_mappability_vars.sh for details
# input:
#     $TN5_KMERS_FILE
#     $N_TN5_PREFERENCE_POSITIONS
# outputs:
#     sorted list of all possible k-mers to $TN5_KMERS_FILE

use strict;
use warnings;

# print all possible kmers of length $nBases in lexicographical order
open my $outH, ">", $ENV{TN5_KMERS_FILE} or die "create_Tn5_kmers error: cannot open $ENV{TN5_KMERS_FILE} for writing\n";
my @bases = ('A', 'C', 'G', 'T');
my $total = 4 ** $ENV{N_TN5_PREFERENCE_POSITIONS};
sub generate_kmers {
    my ($prefix, $length) = @_;
    if ($length == 0) {
        print $outH "$prefix\n";
        return;
    }
    foreach my $base (@bases) {
        generate_kmers($prefix . $base, $length - 1);
    }
}
generate_kmers('', $ENV{N_TN5_PREFERENCE_POSITIONS});
close $outH;
