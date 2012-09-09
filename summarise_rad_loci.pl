#!/usr/bin/env perl

# summarise_RAD_loci.pl
#
# Input: R-style CSV output by RAD_locus_perfect_matches.pl
#        Only works on paired end data
# Output: summary stats for RAD locus coverage
# Author: John Davey john.davey@ed.ac.uk


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

my $tags_filename        = "";
my $fragments_filename   = "";
my $repeats_filename     = "";
my $proportions_filename = "";
my $read_length = 95;
my $options_okay         = GetOptions(
    'tags=s'        => \$tags_filename,
    'frags=s'       => \$fragments_filename,
    'repeats=s'     => \$repeats_filename,
    'proportions=s' => \$proportions_filename,
    'length=i'      => \$read_length,
);

croak
"\nUsage: perl summarise_RAD_loci.pl -t tags_filename -f frags_filename -r repeats_filename -p proportions_filename -l read_length (default 95)\n"
  if !$options_okay;

croak "Please specify a RAD loci CSV file with -t"
  if $tags_filename eq "";

croak "Please specify a fragments file with -f" if $fragments_filename eq "";

croak "Please specify a repeats file with -r" if $repeats_filename eq "";

croak "Please specify a proportions file with -p"
  if $proportions_filename eq "";

#############################################################################
###
### LOAD REPEATS

my %repeats;
if ( $repeats_filename ne "" ) {
    open my $repeats_file, '<', $repeats_filename
      or croak "Can't open repeats file $repeats_filename: $OS_ERROR!\n";

    while ( my $repeats_line = <$repeats_file> ) {
        chomp $repeats_line;
        my ( $chr, $loc, $matches, $reads ) = split /\t/, $repeats_line;
        $repeats{$chr}{$loc}{matches} = $matches;
        $repeats{$chr}{$loc}{reads}   = $reads;
    }

    close $repeats_file;
}

#############################################################################
###
### LOAD FRAGMENTS

open my $fragments_file, '<', $fragments_filename
  or croak "Can't open fragments file $fragments_filename: $OS_ERROR!\n";

my %frag;

my $frag_header = <$fragments_file>;
while ( my $fragment_line = <$fragments_file> ) {
    chomp $fragment_line;
    my ( $chr, $start, $end, $len, $gc, $seq ) = split /\t/, $fragment_line;
    $frag{$chr}{$start}{len}       = $len;
    $frag{$chr}{$start}{genfraggc} = $gc;
    $frag{$chr}{$start}{read1gc}   = gc_content( substr( $seq, 0, $read_length ) );
    $frag{$chr}{$start}{fragseq}   = substr( $seq, 0, 1000 );
    $frag{$chr}{$end}{len}         = $len;
    $frag{$chr}{$end}{genfraggc}   = $gc;
    $frag{$chr}{$end}{read1gc} = gc_content( revcomp( substr( $seq, -$read_length ) ) );
    $frag{$chr}{$end}{fragseq} = revcomp( $seq, -1000 );
}
close $fragments_file;

#############################################################################
###
### LOAD PROPORTIONS

open my $proportions_file, '<', $proportions_filename
  or croak "Can't open proportions file $proportions_filename: $OS_ERROR!\n";

my $prop_header = <$proportions_file>;
my %read_props;
my %frag_props;
while ( my $prop_line = <$proportions_file> ) {
    chomp $prop_line;
    my (
        $sample,       $sample_reads, $sample_reads_prop,
        $sample_frags, $sample_frags_prop
      )
      = split /\t/,
      $prop_line;
    $read_props{$sample} = $sample_reads_prop;
    $frag_props{$sample} = $sample_frags_prop;
}

#############################################################################
###
### LOAD TAGS

open my $tags_file, '<', $tags_filename
  or croak "Can't open tags file $tags_filename: $OS_ERROR!\n";

my %tag;
$tag{loc} = "";

my $header = <$tags_file>;
chomp $header;
my ( $h_chr, $h_loc, $h_taglan, $h_fraglen, $h_gc, $h_totaldepth, @samplenames )
  = split /,/, $header;

my @samplereadprops;
my @samplefragprops;

foreach my $sample (@samplenames) {
    croak "Can't find read proportion for sample $sample\n"
      if ( !defined $read_props{$sample} );
    croak "Can't find frag proportion for sample $sample\n"
      if ( !defined $frag_props{$sample} );
    push @samplereadprops, $read_props{$sample};
    push @samplefragprops, $frag_props{$sample};
}

print
"Chr\tLoc\tGenomeFragLen\tMinTagFragLen\tMaxTagFragLen\tMinTagLenGC\tMeanTagGC\tMaxTagLenGC\tTotalReads\tSample\tPCR\tRep\tSampleFragments\tSampleFragmentsNorm\tSampleReads\tSampleReadsNorm\tMinSampleFragLen\tMeanSampleFragLen\tMaxSampleFragLen\tMinSampleLenGC\tMeanSampleLenGC\tMaxSampleLenGC\n";

while ( my $tags_line = <$tags_file> ) {
    chomp $tags_line;
    my ( $chr, $loc, $taglen, $fraglen, $gc, $totaldepth, @sampledepths ) =
      split /,/, $tags_line;
    next if ( defined $repeats{$chr}{$loc} );

    if ( $loc ne $tag{loc} ) {

        output_tag( \%tag, \@samplenames, \@samplereadprops, \@samplefragprops,
            \%frag, $chr, $loc, $taglen );

        $tag{chr}        = $chr;
        $tag{loc}        = $loc;
        $tag{taglen}     = $taglen;
        $tag{reads}      = 0;
        $tag{minfraglen} = 1000;
        $tag{maxfraglen} = 0;
        $tag{sample}     = ();
        $tag{gc}         = ();

        map { $tag{sample}{$_}{fragments} = 0 } @samplenames;
    }

    $tag{reads} += $totaldepth;
    if ( $gc ne "NA" ) {
        push @{ $tag{gc} }, $gc;
    }

    if ( $totaldepth > 0 ) {
        if ( $fraglen < $tag{minfraglen} ) { $tag{minfraglen} = $fraglen }
        if ( $fraglen > $tag{maxfraglen} ) { $tag{maxfraglen} = $fraglen }
    }

    for my $i ( 0 .. $#sampledepths ) {
        $tag{sample}{ $samplenames[$i] }{depth} += $sampledepths[$i];
        if ( $sampledepths[$i] > 0 ) {
            $tag{sample}{ $samplenames[$i] }{fragments}++;
        }
        for my $j ( 1 .. $sampledepths[$i] ) {
            push @{ $tag{sample}{ $samplenames[$i] }{lengths} }, $fraglen;
        }
    }
}
output_tag( \%tag, \@samplenames, \@samplereadprops, \@samplefragprops, \%frag,
    "", 0, 0 );

close $tags_file;

#############################################################################
###
### SUBROUTINES

sub output_tag {
    my ( $tag_ref, $samplenames_ref, $samplereadprops_ref, $samplefragprops_ref,
        $frag_ref, $chr, $loc, $taglen )
      = @_;

    if ( $tag{loc} ne "" ) {
        if ( $tag_ref->{minfraglen} == 1000 ) { $tag_ref->{minfraglen} = "NA" }
        if ( $tag_ref->{maxfraglen} == 0 )    { $tag_ref->{maxfraglen} = "NA" }

        my $tag_minread2gc = "NA";
        my $tag_maxread2gc = "NA";

        my $fragpos = -1;
        if ( defined $frag_ref->{$chr}{ $tag_ref->{loc} } ) {
            $fragpos = $tag_ref->{loc};
        }
        else {
            my $endpos = $tag_ref->{loc} + $taglen - 1;
            if ( defined $frag_ref->{$chr}{$endpos} ) {
                $fragpos = $endpos;
            }
        }

        if ( $fragpos > -1 ) {
            $tag_minread2gc =
              get_gc_for_len( $frag_ref->{$chr}{$fragpos}{fragseq},
                $tag_ref->{minfraglen} );
            $tag_maxread2gc =
              get_gc_for_len( $frag_ref->{$chr}{$fragpos}{fragseq},
                $tag_ref->{maxfraglen} );
        }

        for my $i ( 0 .. $#{$samplenames_ref} ) {

            # C. elegans library specific, should be ignored for other runs
            my ( $n2, $pcr, $rep ) = split /_/, $samplenames_ref->[$i];
            $pcr = "NA" if (!defined $pcr);
            $rep = "NA" if (!defined $rep);

            my $tag_mean_gc = 0;
            my $sum_gc      = 0;
            if ( defined @{ $tag_ref->{gc} } ) {
                map { $sum_gc += $_ } @{ $tag_ref->{gc} };
                $tag_mean_gc = int( $sum_gc / @{ $tag_ref->{gc} } );
            }
            if ( $tag_mean_gc == 0 ) { $tag_mean_gc = "NA" }

            my $samplereads_norm =
              int( $tag_ref->{sample}{ $samplenames_ref->[$i] }{depth} /
                  $samplereadprops_ref->[$i] );

            my $samplefrags_norm =
              int( $tag_ref->{sample}{ $samplenames_ref->[$i] }{fragments} /
                  $samplefragprops_ref->[$i] );

            my $sample_minfraglen  = 0;
            my $sample_meanfraglen = 0;
            my $sample_maxfraglen  = 0;
            if ( defined @{ $tag{sample}{ $samplenames_ref->[$i] }{lengths} } )
            {
                $sample_minfraglen =
                  min @{ $tag_ref->{sample}{ $samplenames_ref->[$i] }
                      {lengths} };
                $sample_maxfraglen =
                  max @{ $tag_ref->{sample}{ $samplenames_ref->[$i] }
                      {lengths} };
                my $sum_length;
                map { $sum_length += $_ }
                  @{ $tag_ref->{sample}{ $samplenames_ref->[$i] }{lengths} };
                $sample_meanfraglen =
                  int( $sum_length /
                      @{ $tag_ref->{sample}{ $samplenames_ref->[$i] }{lengths} }
                  );
            }

            if ( $sample_minfraglen == 0 )  { $sample_minfraglen  = "NA" }
            if ( $sample_maxfraglen == 0 )  { $sample_maxfraglen  = "NA" }
            if ( $sample_meanfraglen == 0 ) { $sample_meanfraglen = "NA" }

            my $sample_minread2gc  = "NA";
            my $sample_meanread2gc = "NA";
            my $sample_maxread2gc  = "NA";
            if ( $fragpos > -1 ) {
                $sample_minread2gc =
                  get_gc_for_len( $frag_ref->{$chr}{$fragpos}{fragseq},
                    $sample_minfraglen );
                $sample_meanread2gc =
                  get_gc_for_len( $frag_ref->{$chr}{$fragpos}{fragseq},
                    $sample_meanfraglen );
                $sample_maxread2gc =
                  get_gc_for_len( $frag_ref->{$chr}{$fragpos}{fragseq},
                    $sample_maxfraglen );
            }

            print
"$tag_ref->{chr}\t$tag_ref->{loc}\t$tag_ref->{taglen}\t$tag_ref->{minfraglen}\t$tag_ref->{maxfraglen}\t$tag_minread2gc\t$tag_mean_gc\t$tag_maxread2gc\t$tag_ref->{reads}\t$samplenames_ref->[$i]\t$pcr\t$rep\t$tag_ref->{sample}{$samplenames_ref->[$i]}{fragments}\t$samplefrags_norm\t$tag_ref->{sample}{$samplenames_ref->[$i]}{depth}\t$samplereads_norm\t$sample_minfraglen\t$sample_meanfraglen\t$sample_maxfraglen\t$sample_minread2gc\t$sample_meanread2gc\t$sample_maxread2gc\n";
        }
    }

    $tag{chr}        = $chr;
    $tag{loc}        = $loc;
    $tag{taglen}     = $taglen;
    $tag{reads}      = 0;
    $tag{minfraglen} = 1000;
    $tag{maxfraglen} = 0;
    $tag{sample}     = ();
    $tag{gc}         = ();
}

sub get_gc_for_len {
    my ( $seq, $len ) = @_;
    return $len eq "NA" ? "NA" : gc_content( substr( $seq, 0, $len ) );
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
