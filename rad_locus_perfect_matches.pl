#!/usr/bin/env perl

# RAD_locus_perfect_matches.pl
#
# Input: RAD reads in FASTQ format and restriction fragments to map to
# Output: CSV file of aligned RAD loci
#         Read and fragment proportions for each individual sample
#           (paired end data only)
#         Repeat content (paired end data only)
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

use constant QUAL_OFFSET => 33;

# Autoflush output so reporting on progress works
$| = 1;

my $fastq1_filename    = "";
my $fastq2_filename    = "";
my $fragments_filename = "";
my $barcodes_filename  = "";
my $enzyme_w_cutsites  = "C*TGCA*G";
my $read_length        = 100;
my $max_reads          = -1;
my $csv_filename       = "";
my $rpt_filename       = "";
my $prop_filename      = "";

my $options_okay = GetOptions(
    '1fastq=s'        => \$fastq1_filename,
    '2fastq=s'        => \$fastq2_filename,
    'fragments=s'     => \$fragments_filename,
    'length=i'        => \$read_length,
    'barcodes_file=s' => \$barcodes_filename,
    'enzyme=s'        => \$enzyme_w_cutsites,
    'max_reads=i'     => \$max_reads,
    'csv=s'           => \$csv_filename,
    'rpt=s'           => \$rpt_filename,
    'prop=s'          => \$prop_filename,
);

croak
"\nUsage: perl convert_fastq_to_RAD_reads.pl -1 fastq_file_1 -2 fastq_file_2 -f RAD_fragments_file -b barcodes_file -e enzyme -c csv_filename -r repeat_filename -p proportion_filename\n"
  if !$options_okay;

croak "Please specify read 1 file with -1\n" if $fastq1_filename eq "";
croak "Please specify a barcodes file with -b" if $barcodes_filename eq "";
croak "Please specify a file containing RAD fragments with -f"
  if $fragments_filename eq "";

croak "Please specify a filename for CSV output with -c\n" if $csv_filename eq "";
croak "Please specify a filename for repeat output with -r\n" if $rpt_filename eq "";
croak "Please specify a filename for read and fragment proportions output with -p\n" if $prop_filename eq "";

my ( $precut, $centre, $postcut ) = split /\*/, $enzyme_w_cutsites;

croak
"Please specify cut positions on both strands of enzyme with *, eg SbfI=CC*TGCA*GG"
  if ( ( !defined $precut ) || ( !defined $centre ) || ( !defined $postcut ) );

my $enzyme   = $precut . $centre . $postcut;
my $overhang = $centre . $postcut;

my $overhang_length = length $overhang;
my $trim_length     = ( length $enzyme ) - $overhang_length;

#############################################################################
###
### LOAD BARCODES

print STDERR "Loading barcodes...\n";
open my $barcodes_file, '<', $barcodes_filename
  or croak "Can't open barcodes file $barcodes_filename: $OS_ERROR!\n";

my %barcodes;
my %revcomp_barcodes;
my $barcode_length;
while ( my $barcode_line = <$barcodes_file> ) {
    my ( $name, $barcode ) = split /\s+/, $barcode_line;
    $barcodes{$barcode}                    = $name;
    $revcomp_barcodes{ revcomp($barcode) } = $name;
    $barcode_length                        = length $barcode;
}

close $barcodes_file;

my @sample_list = sort values %barcodes;

#############################################################################
###
### LOAD RAD FRAGMENTS

print STDERR "Loading restriction fragments...\n";
open my $fragments_file, '<', $fragments_filename
  or croak "Can't open fragments file $fragments_filename: $OS_ERROR!\n";

my %tags_by_tag;
my %tags_by_chrom;

while ( my $fragments_line = <$fragments_file> ) {
    next if ( $fragments_line =~ /^Chromosome/ );    # Skip header
    chomp $fragments_line;
    my ( $chrom, $start, $end, $frag_length, $gc, $fragment ) = split /\t/,
      $fragments_line;

    my $trim_fragment = substr substr( $fragment, $trim_length ), 0,
      ( length $fragment ) - $trim_length * 2;

    my $frag1 = substr $trim_fragment, 0, 1000;
    my $frag2 =
      length($trim_fragment) > 1000
      ? revcomp( substr $trim_fragment, ( length $trim_fragment ) - 1000 )
      : revcomp($trim_fragment);

    load_fragment( \%tags_by_tag, $frag1, $start, $read_length, $overhang,
        $chrom );
    load_fragment( \%tags_by_tag, $frag2, $end, $read_length, $overhang,
        $chrom );

    if ( $fragment =~ /^$enzyme/ ) {
        $tags_by_chrom{$chrom}{$start}{length} = $frag_length;
        $tags_by_chrom{$chrom}{$start}{seq}    = $fragment;
    }
    if ( $fragment =~ /($enzyme)$/ ) {
        $tags_by_chrom{$chrom}{$end}{length} = $frag_length;
        $tags_by_chrom{$chrom}{$end}{seq}    = revcomp($fragment);
    }

}

close $fragments_file;

#############################################################################
###
### LOAD READS

print STDERR "Loading reads...\n";
open my $fastq1, '<', $fastq1_filename
  or croak "Can't open $fastq1_filename: $OS_ERROR!\n";

my $fastq2;
if ( $fastq2_filename ne "" ) {
    open $fastq2, '<', $fastq2_filename
      or croak "Can't open $fastq2_filename: $OS_ERROR!\n";
}

my $read_count = 0;
my $good_reads = 0;
my $read_found = 1;
my %locations;
my %repeats;

while ($read_found) {
    $read_found = 0;

    my ( $name, $r1, $q1, $r2, $q2 ) = load_read_pair( $fastq1, $fastq2 );

    next if ( !defined $name );

    $read_found = 1;

    $read_count++;
    if ( $read_count % 10000 == 0 ) { print STDERR "." }
    if ( $read_count % 100000 == 0 ) { printf STDERR "%10d", $read_count }
    if ( $read_count % 1000000 == 0 ) {
        printf STDERR ":%10d passed\n", $good_reads;
    }

    my $tag_r1      = substr $r1, $barcode_length;
    my $res_flag    = '-';
    my $sample_bcmm = '-';
    my $matches     = '';

    my $overhang_r1 = substr( $tag_r1, 0, $overhang_length );
    if ( $overhang =~ substr( $tag_r1, 0, $overhang_length ) ) {
        $res_flag = "OK";
    }
    else {
        my $overhang_qual = substr $q1, $barcode_length, $overhang_length;
        $res_flag =
          check_ressite( $r1, $overhang, $centre, $overhang_r1, $overhang_qual,
            $barcode_length, $overhang_length, \%revcomp_barcodes );
    }

    next if ( $res_flag !~ "OK" );

    my $r1_barcode = substr $r1, 0, $barcode_length;
    $sample_bcmm =
      defined( $barcodes{$r1_barcode} ) ? "$barcodes{$r1_barcode}:0" : '-';

    if ( $sample_bcmm eq '-' ) {
        my $barcode_qual = substr $q1, 0, $barcode_length;
        $sample_bcmm =
          find_fuzzy_barcode_match( $r1_barcode, $barcode_qual, \%barcodes );
    }

    next if ( $sample_bcmm eq '-' );

    my ( $samples, $barcode_mismatches ) = split /:/, $sample_bcmm;
    my @samples = split /;/, $samples;
    next if ( @samples > 1 );

    my ( $sample, $substitution ) = split /,/, $samples[0];

    $tag_r1 = substr $tag_r1, $overhang_length;

    my @matches;

    my ( $tag_r1_match, $tag_r1_hd ) =
      defined $tags_by_tag{$tag_r1}
      ? ( $tags_by_tag{$tag_r1}, 0 )
      : undef;

    if ( defined $tag_r1_match ) {
        if ( defined $fastq2 ) {
            my $rev_r2 = revcomp($r2);
            foreach my $uniq_tag ( keys %{$tag_r1_match} ) {
                my $shear_length = 0;

                my @r2_positions;
                while ( $tag_r1_match->{$uniq_tag} =~ /$rev_r2/g ) {
                    push @r2_positions, $+[0];
                }
                my $r2_positions = join ',', @r2_positions;
                $r2_positions .= $r2_positions eq "" ? "-" : "";
                push @matches, "$uniq_tag:$tag_r1_hd:$r2_positions";
            }
        }
        if ( @matches == 0 ) {    # No Read 2 match, but Read 1 still matches
            if ( keys %{$tag_r1_match} == 1 )
            { # If there is only one tag with this Read 1 sequence, assign this read to it
                push @matches,
                  ( keys %{ $tags_by_tag{$tag_r1} } )[0] . ":$tag_r1_hd:-";
            }
        }
    }

    next if ( @matches == 0 );
    if ( @matches == 1 ) {
        my ( $chromosome, $location, $mismatches, $length ) = split /:/,
          $matches[0];
        next if ( ( defined $fastq2 ) && ( $length eq "-" ) );
        $locations{$chromosome}{$location}{$sample}{$length}++;

        $good_reads++;
        last if ( ( $max_reads > 0 ) && ( $good_reads == $max_reads ) );
    }
    else { # Repeat read
        foreach my $match (@matches) {
            my ( $chromosome, $location, $mismatches, $length ) = split /:/,
              $match;
            $repeats{$chromosome}{$location}{reads}++;
            $repeats{$chromosome}{$location}{matches} = scalar @matches;
        }
    }
}
printf STDERR ":%10d passed\n", $good_reads;

close $fastq2 if defined $fastq2;
close $fastq1;

#############################################################################
###
### OUTPUT CSV FILE

print STDERR "Writing CSV file...\n";
open my $csv_file, ">", $csv_filename or croak "Can't open CSV file $csv_filename! $OS_ERROR\n";

if ( $fastq2_filename ne "" ) {
    print $csv_file "Chr,Loc,ResFragLen,ShearFragLen,GC,TotalDepth";
}
else {
    print $csv_file "Chr,Loc,ResFragLen,TotalDepth";
}

map { print $csv_file ",$_" } @sample_list;

my %sample_reads;
my %sample_frags;

print $csv_file "\n";
my $locus_count = 0;
foreach my $chr ( sort keys %tags_by_chrom ) {
    foreach my $loc ( sort { $a <=> $b } keys %{ $tags_by_chrom{$chr} } ) {

        $locus_count++;
        if ( $locus_count % 100 == 0 ) { print STDERR "." }
        if ( $locus_count % 1000 == 0 ) { printf STDERR "%10d", $locus_count }
        if ( $locus_count % 10000 == 0 ) { printf STDERR "\n"; }

        if ( $fastq2_filename ne "" ) {
            foreach my $len ( 100 .. 900 ) {
                my $gc;
                if ( $tags_by_chrom{$chr}{$loc}{length} < $len ) {
                    $gc = 'NA';
                }
                else {
                    my $frag = substr $tags_by_chrom{$chr}{$loc}{seq}, 0, $len;
                    $gc = $frag =~ tr/G//;
                    $gc += $frag =~ tr/C//;
                    $gc = int( ( $gc / $len ) * 100 );
                }

                my $total_depth   = 0;
                my $sample_string = '';
                foreach my $sample (@sample_list) {
                    if ( defined $locations{$chr}{$loc}{$sample}{$len} ) {
                        $sample_string .=
                          ",$locations{$chr}{$loc}{$sample}{$len}";
                        $total_depth += $locations{$chr}{$loc}{$sample}{$len};
                        if ($chr ne "CHROMOSOME_MTDNA") {
                            $sample_reads{$sample}+=$locations{$chr}{$loc}{$sample}{$len};
                            $sample_frags{$sample}++;
                        }
                    }
                    else {
                        $sample_string .= ",0";
                    }
                }
                print $csv_file 
"$chr,$loc,$tags_by_chrom{$chr}{$loc}{length},$len,$gc,$total_depth$sample_string\n";
            }
        }
        else {
            my $total_depth   = 0;
            my $sample_string = '';
            foreach my $sample (@sample_list) {
                if ( defined $locations{$chr}{$loc}{$sample}{'-'} ) {
                    $sample_string .= ",$locations{$chr}{$loc}{$sample}{'-'}";
                    $total_depth += $locations{$chr}{$loc}{$sample}{'-'};
                    if ($chr ne "CHROMOSOME_MTDNA") {
                        $sample_reads{$sample}+=$locations{$chr}{$loc}{$sample}{'-'};
                    }
                }
                else {
                    $sample_string .= ",0";
                }
            }
            print $csv_file 
"$chr,$loc,$tags_by_chrom{$chr}{$loc}{length},$total_depth$sample_string\n";

        }
    }
}

close $csv_file;

print STDERR "\n";

#############################################################################
###
### OUTPUT PROPORTIONS FILE

print STDERR "Writing proportions file...\n";

open my $prop_file, ">", $prop_filename or croak "Can't open proportions file $prop_filename! $OS_ERROR\n";

my $all_read_count = 0;
my $all_frag_count = 0;
foreach my $sample (@sample_list) {
    $all_read_count += $sample_reads{$sample};
    if ($fastq2_filename ne "") {
        $all_frag_count += $sample_frags{$sample};
    }
}

print $prop_file "Sample\tSampleReads\tSampleReadsProp\tSampleFragments\tSampleFragmentsProp\n";
foreach my $sample (@sample_list) {
    my $sample_reads_prop = sprintf "%4.2f", $sample_reads{$sample}/$all_read_count * 100;
    print $prop_file "$sample\t$sample_reads{$sample}\t$sample_reads_prop";

    my $sample_frags_prop = 0;
    if ($fastq2_filename ne "") {
        $all_frag_count += $sample_frags{$sample};
        $sample_frags_prop = sprintf "%4.2f", $sample_frags{$sample}/$all_frag_count * 100;
        print $prop_file "\t$sample_frags{$sample}\t$sample_frags_prop\n";
    }
    else {
        print $prop_file "\t0\t0.00\n";
    }
}

close $prop_file;

#############################################################################
###
### OUTPUT REPEAT FILE

print STDERR "Writing repeats file...\n";

open my $rpt_file, ">", $rpt_filename or croak "Can't open repeat file $rpt_filename! $OS_ERROR\n";

for my $chrom (sort keys %repeats) {
    for my $location (sort {$a<=>$b} keys %{$repeats{$chrom}}) {
        print $rpt_file "$chrom\t$location\t$repeats{$chrom}{$location}{matches}\t$repeats{$chrom}{$location}{reads}\n"
    }
}

close $rpt_file;


print STDERR "Done\n";

#############################################################################
###
### SUBROUTINES

sub load_read_pair {
    my ( $fastq1, $fastq2 ) = @_;

    my $r1name = <$fastq1>;
    my $r1seq  = <$fastq1>;
    my $r1qn   = <$fastq1>;
    my $r1qual = <$fastq1>;

    return if ( ( !$r1name )
        || ( !$r1seq )
        || ( !$r1qn )
        || ( !$r1qual ) );

    chomp $r1name;
    chomp $r1seq;
    chomp $r1qual;

    my ( $r2name, $r2seq, $r2qn, $r2qual );

    if ( defined $fastq2 ) {
        $r2name = <$fastq2>;
        $r2seq  = <$fastq2>;
        $r2qn   = <$fastq2>;
        $r2qual = <$fastq2>;

        return if ( ( !$r2name )
            || ( !$r2seq )
            || ( !$r2qn )
            || ( !$r2qual ) );

        chomp $r2name;
        chomp $r2seq;
        chomp $r2qual;

    }

    if ( $r1name =~ /^@(\S+) (\S+)/xms ) {
        $r1name = $1;
    }

    if ( defined $fastq2 ) {
        if ( $r2name =~ /^@(\S+) (\S+)/xms ) {
            $r2name = $1;
        }

        if ( $r1name ne $r2name ) {
            return $r1name;
        }
    }

    return ( $r1name, $r1seq, $r1qual, $r2seq, $r2qual );

}

sub load_fragment {
    my ( $tag_ref, $frag, $pos, $read_length, $overhang, $chrom ) = @_;
    return if ( $frag !~ /^$overhang/ );
    my $overhang_length = length $overhang;
    my $tag = substr $frag, $overhang_length, $read_length - $overhang_length;

    $tag_ref->{$tag}{"$chrom:$pos"} = $frag;
    return;
}

sub find_fuzzy_barcode_match {
    my ( $r1_barcode, $barcode_qual, $barcodes_ref ) = @_;

    my %hd;
    foreach my $barcode ( sort keys %{$barcodes_ref} ) {
        my $r1hd = hamming( $r1_barcode, $barcode );
        $hd{$r1hd}{$barcode}++;
    }
    my $minhd = min keys %hd;

    if ( keys %{ $hd{$minhd} } == 1 ) {
        my $barcode = ( keys %{ $hd{$minhd} } )[0];
        my $transitions =
          get_transitions( $barcode, $r1_barcode, $barcode_qual );
        return "$barcodes_ref->{$barcode},$transitions:$minhd";
    }
    else {
        my $barcode_list;
        map {
            my $transitions = get_transitions( $_, $r1_barcode, $barcode_qual );
            $barcode_list .= "$barcodes_ref->{$_},$transitions;"
          }
          sort keys %{ $hd{$minhd} };
        chop $barcode_list;
        $barcode_list .= ":$minhd";
        return $barcode_list;
    }

    return '-';
}

sub check_ressite {
    my ( $r1, $overhang, $centre, $overhang_r1, $overhang_qual, $barcode_length,
        $overhang_length, $revcomp_barcodes_ref )
      = @_;

    # Check P1-P2
    if ( $overhang_r1 =~ /^AGATC/ ) {
        return "P2A:-:0";
    }

    # Check P2-P2
    if ( $r1 =~ /^CGTATGCCGT/ ) {
        return "P2P:-:0";
    }

    # Check P1-P1
    if ( $r1 =~ /^[ACGT]{$barcode_length}$centre([ACGT]{$barcode_length})AGA/ )
    {
        if ( defined $revcomp_barcodes_ref->{$1} ) {
            return "P1:-:0";
        }
    }

    if ( $overhang_r1 eq $overhang ) {
        return "OK:-:0";
    }

    my $hd = hamming( $overhang_r1, $overhang );
    if ( $hd <= 2 ) {
        my $transitions =
          get_transitions( $overhang, $overhang_r1, $overhang_qual );
        return "OK:$transitions:$hd";
    }

    for my $i ( 0 .. $barcode_length + $overhang_length ) {
        my $shift_overhang = substr $r1, $i, $overhang_length;
        if ( substr( $r1, $i, $overhang_length ) eq $overhang ) {
            my $shift = $i - $barcode_length;
            return "Sh:-:$shift";
        }
    }

    return '-';
}

sub get_transitions {
    my ( $orig, $seq, $qual ) = @_;

    my @orig = split //, $orig;
    my @seq  = split //, $seq;
    my @qual = split //, $qual;

    my $transitions = "";
    for my $i ( 0 .. $#orig ) {
        if ( $orig[$i] ne $seq[$i] ) {
            my $pos     = $i + 1;
            my $qualnum = ord( $qual[$i] ) - QUAL_OFFSET;
            if ( $qualnum == 2 ) { $qualnum = 0 }
            $transitions .= "$pos($orig[$i]>$seq[$i] $qualnum),";
        }
    }
    chop $transitions;
    return $transitions;
}

sub hamming {
    return ( $_[0] ^ $_[1] ) =~ tr/\001-\255//;
}

sub revcomp {
    my ($seq) = @_;
    $seq = reverse $seq;
    $seq =~ tr/ACGT/TGCA/;
    return $seq;
}
