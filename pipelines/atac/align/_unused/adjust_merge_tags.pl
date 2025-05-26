use strict;
use warnings;

# action:
#     change :0 suffix in QNAME to :1 to reflect merged status
#     drop the rest of the name line added by fastp
# input:
#     merged read FASTQ stream on STDIN
# output: 
#     merged read FASTQ stream on STDOUT with modified QNAME

# working variables
my $nReads = 0;

# run the interleaved pairs
my $lineN = 0;
while(my $line = <STDIN>){
    unless($lineN % 4){
        $nReads++;
        my @name = split(" ", $line); # split drops trailing whitespace when splitting on whitespace, i.e., chomp occurs implicitly
        @name = split(":", $name[0]);
        $name[$#name] = 1;
        $line = join(":", @name)."\n";
    }
    print $line;
    $lineN++;
}

# report read counts
print STDERR "    merged read count = $nReads\n";
