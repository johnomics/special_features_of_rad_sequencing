#!/usr/bin/env perl

# simulate_rad_fragments.pl
#
# Input: reference genome and restriction enzyme
# Output: restriction fragments with start, end, length and GC statistics
# Author: John Davey johnomics@gmail.com


#############################################################################
###                                                                       ###
###                                 CODE                                  ###
###                                                                       ###
#############################################################################

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Data::Dumper;

# Autoflush output so reporting on progress works
$| = 1;


my $in_filename = "";
my $enzyme = "CTGCAG";
my $options_okay = GetOptions(
    'input=s'      => \$in_filename,
    'enzyme=s'     => \$enzyme,
);

croak
"\nUsage: perl simulate_rad_sites.pl -i fastq_file -e enzyme\n"
  if !$options_okay;

croak
"\nPlease specify an input file with -i\nUsage: perl simulate_rad_sites.pl -i fastq_file -e enzyme\n"
  if ($in_filename eq "");

croak
"\nPlease specify an enzyme with -e\nUsage: perl simulate_rad_sites.pl -i fastq_file -e enzyme\n"
  if ($enzyme eq "");


print "Chromosome\tStart\tEnd\tLength\tGC\tSequence\n";

open my $in_file, '<', $in_filename or croak "Can't open input file $in_filename: $OS_ERROR!\n";

my $fragment = "";
my $start_pos = 1;
my $pos = 1;
my $chrom = "";
my $prev_line_end = "";
while (my $in_line = <$in_file>) {
    $in_line = uc $in_line;
    chomp $in_line;
    if ($in_line =~ />(.+)/) {
        if ($fragment ne "") {
            output_fragment($fragment, $prev_line_end, $chrom, $start_pos, $pos);
        }
        $chrom = $1;
        $pos = 1;
        $start_pos = 1;
        $fragment = "";
        $prev_line_end = "";
    }
    else {
        $in_line = $prev_line_end . $in_line;
        while (length $in_line >= length $enzyme) {
            if ($in_line =~ /^$enzyme/i) {
                output_fragment($fragment, $enzyme, $chrom, $start_pos, $pos);
                $in_line = substr $in_line, length $enzyme;
                $start_pos = $pos;
                $pos = $pos + length $enzyme;
                $fragment = $enzyme;
            }
            else {
                $fragment .= substr $in_line, 0, 1;
                $in_line = substr $in_line, 1;
                $pos++;
            }
        }
        $prev_line_end = $in_line;
    }
}

output_fragment($fragment, $prev_line_end, $chrom, $start_pos, $pos);

close $in_file;


sub output_fragment {
    my ($fragment, $end_seq, $chrom, $start_pos, $pos) = @_;
    $fragment .= $end_seq;
    my $end_pos = $pos - 1 + length $end_seq;
    print "$chrom\t$start_pos\t$end_pos\t";
    print length $fragment;
    my $gc = $fragment =~ tr/C//;
    $gc += $fragment =~ tr/G//;
    $gc = int (($gc / length $fragment) * 100);
    print "\t$gc";
    print "\t$fragment\n";
}
