#!/usr/bin/env perl

# infer_rad_loci_statistics.pl
#
# Input: R-style CSV for a single scaffold
# Output: summary statistics for inferred RAD loci
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
use List::Util qw /min max sum/;

# Autoflush output so reporting on progress works
$| = 1;

my $loci_filename = "";
my $ref_filename  = "";
my $tag_length    = 0;
my $output_filename = "";
my $options_okay  = GetOptions(
    'loci=s'       => \$loci_filename,
    'reference=s'  => \$ref_filename,
    'tag_length=i' => \$tag_length,
    'output=s'     => \$output_filename,
);

croak
"\nUsage: perl infer_rad_loci_statistics.pl -l loci_filename -r reference_filename -t tag_length -o output_filename\n"
  if !$options_okay;

croak "Please specify a loci file with -l"
  if $loci_filename eq "";

croak "Please specify a reference with -r" if $ref_filename eq "";

croak
  "Please specify a tag length (ie read length with barcode removed) with -t\n"
  if $tag_length == 0;
  
croak "Please specify an output filename with -o\n" if $output_filename eq "";

#############################################################################
###
### LOAD REFERENCE GENOME

print STDERR "Loading reference $ref_filename...\n";

open my $ref_file, '<', $ref_filename
  or croak "Can't open reference file $ref_filename: $OS_ERROR!\n";

my $seq_name = "";
my %ref;
while ( my $ref_line = <$ref_file> ) {
    chomp $ref_line;
    if ( $ref_line =~ /^>(.+)/ ) {
        $seq_name = ( split " ", $1 )[0];    #Take everything up the first space
    }
    else {
        $ref{$seq_name} .= $ref_line;
    }
}
close $ref_file;

#############################################################################
###
### PROCESS RAD LOCI

print STDERR "Loading loci file $loci_filename...\n";
open my $loci_file, '<', $loci_filename
  or croak "Can't open loci file $loci_filename: $OS_ERROR!\n";

my %tag;

my $header = <$loci_file>;
chomp $header;
my ( $h_chr, $h_loc, $h_dir, $h_seqindex, $h_seq, $h_fraglen, $h_totaldepth,
    @samplenames )
  = split /,/, $header;

while ( my $loci_line = <$loci_file> ) {
    chomp $loci_line;
    my ( $chr, $loc, $dir, $seq_index, $seq, $fraglen, $totaldepth,
        @sampledepths )
      = split /,/, $loci_line;

    next
      if ( $fraglen <= $tag_length )
      ;    # Ignore short fragments likely to contain adapters
    if ( !defined $tag{$chr}{$loc}{$dir}{fragseq} ) {
        if ( $dir eq '+' ) {
            $tag{$chr}{$loc}{$dir}{fragseq} = substr $ref{$chr}, $loc - 1, 1000;
        }
        elsif ( $dir eq '-' ) {
            $tag{$chr}{$loc}{$dir}{fragseq} =
              $loc > 1000
              ? substr( $ref{$chr}, $loc - 1000, 1000 )
              : substr( $ref{$chr}, 0,           $loc );
            $tag{$chr}{$loc}{$dir}{fragseq} =
              revcomp( $tag{$chr}{$loc}{$dir}{fragseq} );
        }
    }

    $tag{$chr}{$loc}{$dir}{allele}{$seq_index}{seq} = $seq;
    $tag{$chr}{$loc}{$dir}{allele}{$seq_index}{frags}{$fraglen}{reads} =
      $totaldepth;
    $tag{$chr}{$loc}{$dir}{totaldepth} += $totaldepth;
    for my $i ( 0 .. $#samplenames ) {
        $tag{$chr}{$loc}{$dir}{allele}{$seq_index}{frags}{$fraglen}{samples}
          { $samplenames[$i] } = $sampledepths[$i];
        if ( $sampledepths[$i] > 0 ) {
            $tag{$chr}{$loc}{$dir}{resfraglen}{ $samplenames[$i] } = 0;
            $tag{$chr}{$loc}{$dir}{allele}{$seq_index}{sample}
              { $samplenames[$i] }{reads} += $sampledepths[$i];
            $tag{$chr}{$loc}{$dir}{allele}{$seq_index}{sample}
              { $samplenames[$i] }{frags}++;
            $tag{$chr}{$loc}{$dir}{allele}{$seq_index}{reads} +=
              $sampledepths[$i];
        }
    }
}

close $loci_file;

#############################################################################
###
### CALCULATE RESTRICTION FRAGMENT LENGTHS AND LOCUS DEPTHS

print STDERR
  "Calculating sample restriction fragment lengths and locus depths...\n";

foreach my $chr ( sort keys %tag ) {

    my @tagloc = sort { $a <=> $b } keys %{ $tag{$chr} };

    for my $i ( 0 .. $#tagloc ) {
        foreach my $dir ( keys %{ $tag{$chr}{ $tagloc[$i] } } ) {
            foreach my $sample (
                sort keys %{ $tag{$chr}{ $tagloc[$i] }{$dir}{resfraglen} } )
            {
                my $match_dir  = $dir eq '+' ? '-' : '+';
                my $search_inc = $dir eq '+' ? 1   : -1;

                my $candidate_loc = $i + $search_inc;
                my $found_match   = 0;
                while (( $candidate_loc >= 0 )
                    && ( $candidate_loc < $#tagloc ) )
                {
                    if (
                        defined $tag{$chr}{ $tagloc[$candidate_loc] }
                        {$match_dir} )
                    {
                        if ( $tagloc[$candidate_loc] ==
                            $tagloc[$i] + 5 * $search_inc )
                        {
                            $candidate_loc += $search_inc;
                            next;
                        }
                        if (
                            defined $tag{$chr}{ $tagloc[$candidate_loc] }
                            {$match_dir}{resfraglen}{$sample} )
                        {
                            $tag{$chr}{ $tagloc[$i] }{$dir}{resfraglen}
                              {$sample} =
                              ( $dir eq '+' )
                              ? $tagloc[$candidate_loc] - $tagloc[$i] + 1
                              : $tagloc[$i] - $tagloc[$candidate_loc] + 1;
                            $tag{$chr}{ $tagloc[$i] }{$dir}{neighbour}
                              {$sample} = $tagloc[$candidate_loc];
                            last;
                        }
                    }
                    $candidate_loc += $search_inc;
                }

                foreach my $allele (
                    keys %{ $tag{$chr}{ $tagloc[$i] }{$dir}{allele} } )
                {
                    if (
                        defined $tag{$chr}{ $tagloc[$i] }{$dir}{allele}{$allele}
                        {sample}{$sample} )
                    {
                        $tag{$chr}{ $tagloc[$i] }{$dir}{locusreads}{$sample} +=
                          $tag{$chr}{ $tagloc[$i] }{$dir}{allele}{$allele}
                          {sample}{$sample}{reads};
                        $tag{$chr}{ $tagloc[$i] }{$dir}{locusfrags}{$sample} +=
                          $tag{$chr}{ $tagloc[$i] }{$dir}{allele}{$allele}
                          {sample}{$sample}{frags};
                    }
                }
            }
        }
    }
}

#############################################################################
###
### OUTPUT RAD LOCUS STATISTICS

print STDERR "Writing output...\n";

open my $out_file, '>', $output_filename or croak "Can't open $output_filename: $OS_ERROR\n";

print $out_file 
"Chr\tLoc\tDir\tAllele\tSeq\tMinTagFragLen\tMaxTagFragLen\tMinTagLenGC\tMaxTagLenGC\tTotalLocusReads\tTotalAlleleReads\tSample\tSampleLocusReads\tSampleLocusFragments\tSampleAlleleReads\tSampleAlleleFragments\tSampleAlleleReadsProp\tSampleAlleleFragsProp\tSampleResFragLen\tSampleMatchingTag\tMinSampleFragLen\tMeanSampleFragLen\tMaxSampleFragLen\tMinSampleLenGC\tMeanSampleLenGC\tMaxSampleLenGC\n";

foreach my $chr ( sort keys %tag ) {
    foreach my $loc ( sort { $a <=> $b } keys %{ $tag{$chr} } ) {
        foreach my $dir ( keys %{ $tag{$chr}{$loc} } ) {
            my $fragseq = $tag{$chr}{$loc}{$dir}{fragseq};
            foreach my $allele ( sort { $a <=> $b }
                keys %{ $tag{$chr}{$loc}{$dir}{allele} } )
            {

                my $tag_ref = $tag{$chr}{$loc}{$dir}{allele}{$allele};

                my $mintagfraglen = min( keys %{ $tag_ref->{frags} } );
                my $maxtagfraglen = max( keys %{ $tag_ref->{frags} } );
                my $mintagfraglengc =
                  get_gc_for_len( $fragseq, $mintagfraglen );
                my $maxtagfraglengc =
                  get_gc_for_len( $fragseq, $maxtagfraglen );

                foreach my $sample (@samplenames) {
                    print $out_file "$chr\t$loc\t$dir\t$allele\t$tag_ref->{seq}";
                    print $out_file 
"\t$mintagfraglen\t$maxtagfraglen\t$mintagfraglengc\t$maxtagfraglengc";
                    print $out_file 
"\t$tag{$chr}{$loc}{$dir}{totaldepth}\t$tag{$chr}{$loc}{$dir}{allele}{$allele}{reads}";

                    if ( !defined $tag{$chr}{$loc}{$dir}{locusreads}{$sample} )
                    {
                        $tag{$chr}{$loc}{$dir}{locusreads}{$sample} = "NA";
                    }

                    if ( !defined $tag{$chr}{$loc}{$dir}{locusfrags}{$sample} )
                    {
                        $tag{$chr}{$loc}{$dir}{locusfrags}{$sample} = "NA";
                    }

                    print $out_file  "\t$sample";
                    print $out_file 
"\t$tag{$chr}{$loc}{$dir}{locusreads}{$sample}\t$tag{$chr}{$loc}{$dir}{locusfrags}{$sample}";

                    if ( defined $tag_ref->{sample}{$sample}{reads} ) {
                        if ( !defined $tag{$chr}{$loc}{$dir}{neighbour}
                            {$sample} )
                        {
                            $tag{$chr}{$loc}{$dir}{neighbour}{$sample} = "NA";
                        }

                        if (
                            $tag{$chr}{$loc}{$dir}{resfraglen}{$sample} eq "0" )
                        {
                            $tag{$chr}{$loc}{$dir}{resfraglen}{$sample} = "NA";
                        }

                        print $out_file 
"\t$tag_ref->{sample}{$sample}{reads}\t$tag_ref->{sample}{$sample}{frags}";

                        my $sample_allele_reads_prop = sprintf "%3.2f",
                          $tag_ref->{sample}{$sample}{reads} /
                          $tag{$chr}{$loc}{$dir}{locusreads}{$sample};
                        my $sample_allele_frags_prop = sprintf "%3.2f",
                          $tag_ref->{sample}{$sample}{frags} /
                          $tag{$chr}{$loc}{$dir}{locusfrags}{$sample};
                        print $out_file 
"\t$sample_allele_reads_prop\t$sample_allele_frags_prop";

                        print $out_file 
"\t$tag{$chr}{$loc}{$dir}{resfraglen}{$sample}\t$tag{$chr}{$loc}{$dir}{neighbour}{$sample}";

                        my @sample_fraglens;
                        foreach my $fraglen ( keys %{ $tag_ref->{frags} } ) {
                            if (
                                (
                                    defined $tag_ref->{frags}{$fraglen}{samples}
                                    {$sample}
                                )
                                && ( $tag_ref->{frags}{$fraglen}{samples}
                                    {$sample} > 0 )
                              )
                            {
                                push @sample_fraglens, $fraglen;
                            }
                        }
                        my $minsamplefraglen = min(@sample_fraglens);
                        my $meansamplefraglen =
                          @sample_fraglens > 0
                          ? int( sum(@sample_fraglens) / @sample_fraglens )
                          : 0;
                        my $maxsamplefraglen = max(@sample_fraglens);

                        my $minsamplefraggc =
                          get_gc_for_len( $fragseq, $minsamplefraglen );
                        my $meansamplefraggc =
                          get_gc_for_len( $fragseq, $meansamplefraglen );
                        my $maxsamplefraggc =
                          get_gc_for_len( $fragseq, $maxsamplefraglen );

                        print $out_file 
"\t$minsamplefraglen\t$meansamplefraglen\t$maxsamplefraglen";
                        print $out_file 
"\t$minsamplefraggc\t$meansamplefraggc\t$maxsamplefraggc";
                    }
                    else {
                        print $out_file 
                          "\t0\t0\t0.00\t0.00\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
                    }
                    print $out_file "\n";
                }
            }
        }
    }
}

close $out_file;

print STDERR "Done\n";

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
