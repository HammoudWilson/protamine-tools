#!/usr/bin/perl
use strict;
use warnings;

# this is a streamlined version of the segment utility found here:
# https://git.umms.med.umich.edu/wilson_lab_public/utilities/-/blob/master/segment/segment

# Variable designations in this script follow as closely as possible 
# the description of the Viterbi algorithm from:
#    http://www.cs.brown.edu/research/ai/dynamics/tutorial/Documents/HiddenMarkovModels.html
#
# 'segment' solves the Viterbi algorithm, i.e. determines the most likely sequence
# of N states over T observations, for a provided Hidden Mardov Model that obeys
# the following properties:
#     state values are indexed in an integer series from 0 to (N-1)
#     observation values are indexed in an integer series from 0 to (M-1)
#
# Usage:
#     perl segment.pl
# where the following environment variable must be set:
#     PERSISTENCE = P(state_t+1 = state_t), i.e. probability of remaining in state
#
# Data are input and output in a stream as follows:
#     stdin  = emissionProbs [optional sequence instance identifiers ...]
#     stdout = inferredState [optional sequence instance identifiers ...]
# where:
#     emissionProbs is a comma-delimited list of emission probabilities in format:
#         P(instance | state=0),P(instance | state=1)[...,P(instance | state=N-1)]
#     inferredState is a value in state set 0 to (N-1)
# and where:
#     all input probabilities, including PERSISTENCE and emissionProbs, are provided as log values
#
# e.g., providing a data stream on stdin such as:
#     0.12,0.88 chr1 1250 +
#     0.35,0.65 chr1 1300 +
#     ...
# might result in a data stream on stdout such as:
#     1 chr1 1250 +
#     1 chr1 1300 +
#     ...

# initialize the data structures
my (@A, @B, @pi, @prevDelta, @phi, @segIDs, @sigma, @i_star);

# collect the bins and observations
while(my $bin = <STDIN>){
    $bin =~ m/^(\S+)(.*)/ or die "segment.pl error: format error in sequence line:\n\t$bin\n";
    my ($obs, $segID) = ($1, $2);
    $obs = [split(",", $obs)];
    defined $segID or $segID = "";
    push @sigma,  $obs;   #sequence of actual observations OR bin emission probabilities 
    push @segIDs, $segID; #undefined list of sequence identifiers supplied by caller
}

# determine state and observation counts
my $N = @{$sigma[0]};  #number of states
$N >= 2 or die "segment.pl error: fewer than two states provided\n"; 
my $N_ = $N - 1; #max 0-referenced state index
my $N__ = $N_ - 1; 
my $T = @sigma; #length of the sequence of observations
my $T_ = $T - 1; #max 0-referenced observation index

# start probabilities derived from ZERO_PROB
my $ZERO_PROB = 0.5;
my $non_zero_prob = (1 - $ZERO_PROB) / $N_; # split non-zero starting probability equally among all non-zero states
foreach my $i(0..$N_){ $pi[$i] = log($i == 0 ? $ZERO_PROB : $non_zero_prob) }

# transition probabilities 
my $PERSISTENCE = exp($ENV{PERSISTENCE}); # since expected to be set as a log value
my $trans_zero_nonZero = (1 - $PERSISTENCE) / $N_;
my $trans_nonZero_zero = (1 - $PERSISTENCE) * ($ZERO_PROB / ($ZERO_PROB + ($N__ * $non_zero_prob))); 
my $trans_nonZero_nonZero = 1 - $PERSISTENCE - $trans_nonZero_zero;
my $out_of_state = (1 - $PERSISTENCE) / $N_;
foreach my $i(0..$N_){
    foreach my $j(0..$N_){
        my $p = 0;
        if($i == $j){
            $p = $PERSISTENCE;  #probability of remaining in state
        } else {
            $p = $out_of_state; #probability of changing state, unweighted
        }  
        $A[$i][$j] = log($p);
    }
}

# 1. initialization
my $logZero = -1E100;
foreach my $i(0..$N_){ 
    push @prevDelta, $pi[$i] + $sigma[0][$i]; #probs are added since using log probs
    push @{$phi[0]}, -9999; #placeholder, value never used  
}

# 2. recursion
foreach my $t(1..$T_){
    my ($o_t, @delta) = ($sigma[$t]);
    foreach my $j(0..$N_){ 
        my ($max_i, $argmax_i) = $logZero;
        foreach my $i(0..$N_){ 
            my $delta = $prevDelta[$i] + $A[$i][$j] + $$o_t[$j];
            $max_i < $delta and $max_i = $delta and $argmax_i = $i;
        }
        push @delta, $max_i;
        push @{$phi[$t]}, $argmax_i; 
    }
    @prevDelta = @delta;
}

# 3. termination
my ($max_i, $argmax_i) = ($logZero);
foreach my $i(0..$N_){ 
    my $delta = $prevDelta[$i];
    $max_i < $delta and $max_i = $delta and $argmax_i = $i;
}
my $p_star = $max_i;
# print STDERR "Viterbi path log probability = $p_star\n";
unshift @i_star, $argmax_i;

# 4. reconstruction
foreach my $t(reverse(1..$T_)){
    unshift @i_star, $phi[$t][$i_star[0]];
}

# print the results
foreach my $t(0..$T_){
    print "$i_star[$t]$segIDs[$t]\n";
}
