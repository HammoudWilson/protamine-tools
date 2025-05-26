use strict;
use warnings;

# action:
#     trim SEQ and QUAL to 150bp, i.e., drop the inconsistently-present 151st base
#     trim trailing data from QNAME (not critical but done)
# input:
#     FASTQ stream from one read on STDIN
# output: 
#     FASTQ stream from one read on STDOUT with trimmed SEQ and QUAL

# constants
use constant READ_LEN => 150;

# working variables
my ($readN) = @ARGV;
my $nReads = 0;

# run the reads
while (my $name = <STDIN>){
    my $seq     = <STDIN>;
    my $discard = <STDIN>;
    my $qual    = <STDIN>;
    chomp $name;
    chomp $seq;
    chomp $qual;
    my @name = split(" ", $name); 
    print join("\n", 
        $name[0],
        substr($seq, 0, READ_LEN), 
        "+", 
        substr($qual, 0, READ_LEN)
    ), "\n";
    $nReads++;
}

# report read counts
print STDERR "    read$readN count = $nReads\n";
