#!/usr/bin/env perl

# infer_rad_loci_per_scaffold_from_bam.pl
#
# Input: one BAM file covering multiple individuals
# Output: RAD read pileup stats in CSV format, one file for each scaffold
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
use List::Util qw /min max/;

# Autoflush output so reporting on progress works
$| = 1;

my $bam_filename   = "";
my $output_folder  = "";
my $enzyme_overhang = "";

my $options_okay =
  GetOptions( 'bam=s' => \$bam_filename, 'enzyme_overhang=s' => \$enzyme_overhang, 'output=s' => \$output_folder);

croak
"\nUsage: perl infer_rad_loci_per_scaffold_from_bam.pl -b bam_file -e enzyme_overhang -o output_folder\n"
  if !$options_okay;

croak "Please specify a BAM file with -b\n"   if ( $bam_filename   eq "" );
croak "Please specify a folder for output with -o\n" if ($output_folder eq "");
croak "Please specify a restriction enzyme overhang with -e (eg for PstI, cut site CTGCAG, overhang is TGCAG)\n" if ($enzyme_overhang eq "");

my %frag;
my %sample;

my %valid_flag = ( 147 => 0, 99 => 0, 163 => 0, 83 => 0 );
my $cur_scf = "";

#############################################################################
###
### PROCESS BAM FILE

print STDERR "Loading BAM file...\n";

open my $bam_file, "-|", "samtools view $bam_filename" or croak "Can't open $bam_filename for viewing:$OS_ERROR\n";

while ( my $bam_line = <$bam_file> ) {
    chomp $bam_line;
    my (
        $read_name, $flag,     $scf,     $pos, $mapq, $cigar,
        $pair_site, $pair_loc, $fraglen, $seq, $qual, @f
    ) = split /\t/, $bam_line;
    if ($scf ne $cur_scf) {
        if ($cur_scf ne "") {
            if (keys %frag > 0) {
                output_scaffold($cur_scf, $output_folder, $bam_filename, \%frag, \%sample);
            }
	    undef %frag;
            undef %sample;
        }
        $cur_scf = $scf;
    }

    next if ( $cigar eq "*" );

    next if ( !defined( $valid_flag{$flag} ) );

    my $read = ( $flag & 64 ) ? 1 : 2;

    # Skip paired end reads (can get everything from read 1 fraglen)
    next if ($read == 2);

#    Discard short fragments less than read length
#    next if ( abs($fraglen) < length($seq) );

    my $sample_name   = "NA";

    for my $field_group (@f) {
        my ( $field, $type, $value ) = split /:/, $field_group;

        if ( $field eq "RG" ) {
            $sample_name = $value;

            # Details for Heliconius BAM, shouldn't affect other BAM files
            if ( $sample_name =~ /PstI\.(.+)\.110802/ ) {
                $sample_name = $1;

                if ( $sample_name eq "F1F" ) { $sample_name = "F1Mo" }
                if ( $sample_name eq "F1M" ) { $sample_name = "F1Fa" }
                if ( $sample_name eq "F0R" ) { $sample_name = "F0GM" }
            }
        }
    }

    if ( $sample_name ne "NA" ) {
        $sample{$sample_name}++;
    }

    $fraglen=abs($fraglen);
    my $frag_dir     = ' ';
    my $frag_start = $pos + cigar_offset_start($cigar);

    my $revcomp = $flag & 16;
    if ( $revcomp ) {
        $frag_start +=
          length($seq) + cigar_offset_end($cigar);
        $seq = revcomp($seq);
        $frag_dir = '-';
    }
    else {
        $frag_start -= 1;
        $frag_dir = '+';
    }
    
    if ( ( $frag_dir ne ' ' ) && ( $seq =~ /^$enzyme_overhang/ ) ) {
        $frag{$frag_start}{$frag_dir}{readseqs}{$seq}
          {totalreads}++;
        $frag{$frag_start}{$frag_dir}{readseqs}{$seq}{samples}
          {$sample_name}{$fraglen}++;
    }
}
output_scaffold($cur_scf, $output_folder, $bam_filename, \%frag, \%sample);

print STDERR "\n";
close $bam_file;

#############################################################################
###
### SUBROUTINES

sub output_scaffold {
    my ($scf, $output_folder, $bam_filename, $frag_ref, $sample_ref) = @_;

    my $header =  "Scf,Loc,Dir,SeqIndex,Seq,ReadFragLen,TotalDepth";
    map { $header .= ",$_" } sort keys %{$sample_ref};
    $header .= "\n";
    my $file_begun = 0;
    my $scaffold_file;
    print STDERR "Writing scaffold $scf...\n";
    foreach my $pos ( sort { $a <=> $b } keys %{ $frag_ref } ) {

        foreach my $dir ( sort keys %{ $frag_ref->{$pos} } ) {
            my $seq_index = 0;

            foreach my $seq (
                sort {
                    $frag_ref->{$pos}{$dir}{readseqs}{$b}
                      {totalreads} <=> $frag_ref->{$pos}{$dir}{readseqs}{$a}
                      {totalreads}
                } keys %{ $frag_ref->{$pos}{$dir}{readseqs} }
              )
            {
                next
                  if (
                    keys %{ $frag_ref->{$pos}{$dir}{readseqs}{$seq}{samples} }
                    == 1 );
                $seq_index++;
                foreach my $fraglen ( 25 .. 900 ) {
                    my $output = "$scf,$pos,$dir,$seq_index,$seq,$fraglen";
                    my $sample_string = "";
                    my $total_reads   = 0;
                    foreach my $sample ( sort keys %{$sample_ref} ) {
                        if (
                            defined $frag_ref->{$pos}{$dir}{readseqs}{$seq}
                            {samples}{$sample}{$fraglen} )
                        {
                            $sample_string .=
    ",$frag_ref->{$pos}{$dir}{readseqs}{$seq}{samples}{$sample}{$fraglen}";
                            $total_reads +=
                              $frag_ref->{$pos}{$dir}{readseqs}{$seq}{samples}
                              {$sample}{$fraglen};
                        }
                        else {
                            $sample_string .= ",0";
                        }
                    }
                    next if ( $total_reads == 0 );
                    if (!$file_begun) {
                        open $scaffold_file, ">", "$output_folder/$bam_filename\.$scf\.csv" or croak "Can't open output file for scaffold $scf: $OS_ERROR!\n";
                        print $scaffold_file $header;
                        $file_begun++;
                    }
                    print $scaffold_file $output . ",$total_reads$sample_string\n";
                }
            }
        } 
    }
    if ($file_begun) {close $scaffold_file;}
}


sub cigar_offset_start {
    my ($cigar) = @_;
    my $offset = 0;

    if ( $cigar =~ /^([0-9]+)I/ ) {
        return -$1;
    }
    return 0;
}

sub cigar_offset_end {
    my ($cigar) = @_;
    my $offset = 0;

    return 0 if ( $cigar eq "*" );

    while ( $cigar =~ s/([0-9]+)([MIDNSHPX=])// ) {
        my $len  = $1;
        my $type = $2;
        if ( $type !~ /[MID]/ ) {
            print STDERR "Warning: CIGAR type $type ignored\n";
            next;
        }

        # Subtract insertions from reference end position but add deletions
        $offset += $type eq "I" ? -$len : $type eq "D" ? $len : 0;
    }
    return $offset;
}

sub gc_content {
    my ($fragment) = @_;
    my $gc = $fragment =~ tr/C//;
    $gc += $fragment =~ tr/G//;
    $gc = int( ( $gc / length $fragment ) * 100 );
    return $gc;
}

sub revcomp {
    my ($seq) = @_;
    $seq = reverse $seq;
    $seq =~ tr/ACGT/TGCA/;
    return $seq;
}
