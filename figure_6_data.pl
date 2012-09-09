#!/usr/bin/env perl

# figure_6_data.pl
#
# Input: Heliconius TSV files with haplotypes
# Output: Missing genotype data for Figure 6
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

my @genotypes = ( "AB", "AC", "CB", "CC" );
my %parents = (
    "F1Mo" => "AC",
    "F1Fa" => "BC",
    "F0GM" => "AB",
);

my %hapdir = (
    AB => {A=>-1,B=>1},
    AC => {A=>-1,C=>1},
    CB => {C=>-1,B=>1},
    BC => {B=>-1,C=>1},
    CC => {C=>-1},
);

my $tsv_filename = "";

my $options_okay = GetOptions(
    'tsv=s'        => \$tsv_filename,
);

croak
"\nUsage: perl rad_missing_data_plot.pl -t tsv_file\n"
  if !$options_okay;

croak "Please specify a TSV file with -t\n" if ( $tsv_filename eq "" );


#############################################################################
###
### LOAD TSV FILE WITH HAPLOTYPES

print STDERR "Loading TSV file...\n";

my %rad_loci;

open my $tsv_file, '<', $tsv_filename
  or croak "Can't open $tsv_filename: $OS_ERROR!\n";
my $tsv_header = <$tsv_file>;
chomp $tsv_header;
my @tsv_fields = split /\t/, $tsv_header;
my %tsv_field;
for my $i ( 0 .. $#tsv_fields ) {
    $tsv_field{$i} = $tsv_fields[$i];
}

my $genotype_index = 0;
my $prev_sample    = 0;
my $first_allele   = 0;
my $cur_allele     = "";
my $cur_site       = "";

while ( my $tsv_line = <$tsv_file> ) {
    chomp $tsv_line;
    my @f = split /\t/, $tsv_line;
    my %f;
    for my $i ( 0 .. $#f ) {
        $f{ $tsv_field{$i} } = $f[$i];
    }

    next if (($f{MoFragLen} eq "NA") && ($f{FaFragLen} eq "NA") && ($f{SharedFragLen} eq "NA")); # Ignore samples with no fragment length information
    next if $f{Allele} > 3;    #Ignore alleles 4 and above

    next if ($f{Sample} eq "F0GM");
    my $locid = "$f{Chr}.$f{Loc}.$f{Dir}";

    # Get parental genotype
    if ( $f{Sample} =~ /^F/ ) {
        $rad_loci{$locid}{sample}{ $f{Sample} }{gt} =
          $parents{ $f{Sample} };
    }
    else {

        # Get offspring genotype
        if ( $prev_sample > $f{Sample} ) {
            $genotype_index++;
        }
        $prev_sample = $f{Sample};

        my $full_allele = "$f{Chr}$f{Loc}$f{Dir}$f{Allele}";
        if ( $cur_allele ne $full_allele ) {
            $genotype_index = 0;
            $cur_allele     = $full_allele;
        }
        my $full_site = "$f{Chr}$f{Loc}$f{Dir}";
        if ( $cur_site ne $full_site ) {
            $cur_site     = $full_site;
            $first_allele = $f{Allele};
        }
        $rad_loci{$locid}{sample}{ $f{Sample} }{gt} =
          $genotypes[$genotype_index];
    }

    if (defined $rad_loci{$locid}{fraglen} && $rad_loci{$locid}{fraglen} ne $f{SharedFragLen}) {
        delete $rad_loci{$locid};
        next;
    }

    $rad_loci{$locid}{dir} = $f{Dir};
    $rad_loci{$locid}{fraglen} = $f{SharedFragLen};
    foreach my $hap ($f{MoHap}, $f{FaHap}) {
        next if ($hap eq "NA");
        $rad_loci{$locid}{haplotypes}{$hap} = $f{Allele};
        my $short_hap = substr $hap, -1, 1;
        $rad_loci{$locid}{shorthap}{$short_hap} = $f{Allele};
    }

    if ($f{AlleleCount} == 2) {
        $rad_loci{$locid}{sample}{$f{Sample}}{1}{coverage} = $f{SampleAlleleReads}/2;
        $rad_loci{$locid}{sample}{$f{Sample}}{2}{coverage} = $f{SampleAlleleReads}/2;
        $rad_loci{$locid}{sample}{$f{Sample}}{1}{seqid} = $f{Allele};
        $rad_loci{$locid}{sample}{$f{Sample}}{2}{seqid} = $f{Allele};
        $rad_loci{$locid}{sample}{$f{Sample}}{1}{hap} = substr $f{MoHap}, -1, 1;
        $rad_loci{$locid}{sample}{$f{Sample}}{2}{hap} = substr $f{FaHap}, -1, 1;
        $rad_loci{$locid}{sample}{$f{Sample}}{1}{hapdir} = $hapdir{$rad_loci{$locid}{sample}{$f{Sample}}{gt}}{$rad_loci{$locid}{sample}{$f{Sample}}{1}{hap}};
        $rad_loci{$locid}{sample}{$f{Sample}}{2}{hapdir} = $hapdir{$rad_loci{$locid}{sample}{$f{Sample}}{gt}}{$rad_loci{$locid}{sample}{$f{Sample}}{2}{hap}};
        if ($rad_loci{$locid}{sample}{$f{Sample}}{1}{hap} eq "C" && $rad_loci{$locid}{sample}{$f{Sample}}{2}{hap} eq "C") {
            $rad_loci{$locid}{sample}{$f{Sample}}{1}{hapdir} = "-1";
            $rad_loci{$locid}{sample}{$f{Sample}}{2}{hapdir} = "1";
        }
    }
    elsif ($f{AlleleCount} == 1) {
        my $allele = defined $rad_loci{$locid}{sample}{$f{Sample}}{1} ? 2 : 1;
        $rad_loci{$locid}{sample}{$f{Sample}}{$allele}{coverage} = $f{SampleAlleleReads};
        $rad_loci{$locid}{sample}{$f{Sample}}{$allele}{seqid} = $f{Allele};
        my $hap = $f{MoHap} ne "NA" ? $f{MoHap} : $f{FaHap};
        $rad_loci{$locid}{sample}{$f{Sample}}{$allele}{hap} = substr $hap, -1, 1;
        $rad_loci{$locid}{sample}{$f{Sample}}{$allele}{hapdir} = $hapdir{$rad_loci{$locid}{sample}{$f{Sample}}{gt}}{$rad_loci{$locid}{sample}{$f{Sample}}{$allele}{hap}};
    }
}
close $tsv_file;

#############################################################################
###
### DEFINE COLOURS FOR ALLELES

my %fraglen;
foreach my $locid (keys %rad_loci) {
    if (keys %{$rad_loci{$locid}{haplotypes}} < 4) {
        delete $rad_loci{$locid};
        next;
    }
    $fraglen{$rad_loci{$locid}{fraglen}}{$locid}++;
    
    my $hapcolour = 1;
    foreach my $hap ("C", "A","B") {
        $hapcolour++ if ($rad_loci{$locid}{shorthap}{$hap} ne $rad_loci{$locid}{shorthap}{C});
        $rad_loci{$locid}{hapcolour}{$hap} = $hapcolour;
    }
}

#############################################################################
###
### COUNT INDIVIDUALS WITH MISSING ALLELES FOR EACH LOCUS

foreach my $locid (sort {$rad_loci{$a}{fraglen}<=>$rad_loci{$b}{fraglen}} keys %rad_loci) {
    $rad_loci{$locid}{missing} = 0;
    $rad_loci{$locid}{hethom}  = 0;
    foreach my $sample (sort {$rad_loci{$locid}{sample}{$a}{gt} cmp $rad_loci{$locid}{sample}{$b}{gt} || (($a =~ /^F/ || $b =~ /^F/) ? $a cmp  $b : $a <=> $b)} keys %{$rad_loci{$locid}{sample}}) {
        if (defined $rad_loci{$locid}{sample}{$sample}{2}) {
            if (($rad_loci{$locid}{sample}{$sample}{1}{coverage} == 0) && ($rad_loci{$locid}{sample}{$sample}{2}{coverage} == 0)) {
                $rad_loci{$locid}{missing}++;
            }
            elsif (($rad_loci{$locid}{sample}{$sample}{1}{coverage} == 0) || ($rad_loci{$locid}{sample}{$sample}{2}{coverage} == 0)) {
                $rad_loci{$locid}{hethom}++;
            }
        }
        else {
            $rad_loci{$locid}{missing}++ if $rad_loci{$locid}{sample}{$sample}{1}{coverage} == 0;
        }
    }
}

#############################################################################
###
### OUTPUT LOCI MISSING DATA SUMMARY

print "Locus\tFragLen\tDir\tSample\tGenotype\tAllele\tHaplotype\tHapDir\tHapCol\tSeqId\tDepth\tMissing\tHetHom\n";
foreach my $locid (sort {$rad_loci{$a}{fraglen}<=>$rad_loci{$b}{fraglen}} keys %rad_loci) {
    foreach my $sample (sort {$rad_loci{$locid}{sample}{$a}{gt} cmp $rad_loci{$locid}{sample}{$b}{gt} || (($a =~ /^F/ || $b =~ /^F/) ? $a cmp  $b : $a <=> $b)} keys %{$rad_loci{$locid}{sample}}) {
        foreach my $allele (1, 2) {
            print "$locid\t$rad_loci{$locid}{fraglen}\t$rad_loci{$locid}{dir}\t$sample\t$rad_loci{$locid}{sample}{$sample}{gt}\t$allele\t$rad_loci{$locid}{sample}{$sample}{$allele}{hap}\t$rad_loci{$locid}{sample}{$sample}{$allele}{hapdir}\t$rad_loci{$locid}{hapcolour}{$rad_loci{$locid}{sample}{$sample}{$allele}{hap}}\t$rad_loci{$locid}{sample}{$sample}{$allele}{seqid}\t$rad_loci{$locid}{sample}{$sample}{$allele}{coverage}\t$rad_loci{$locid}{missing}\t$rad_loci{$locid}{hethom}\n";
        }
    }
}