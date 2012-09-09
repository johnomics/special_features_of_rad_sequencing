#!/usr/bin/env perl

# infer_rad_loci_haplotypes.pl
#
# Input: unique tags TSV file
# Output: haplotype occurrences across one scaffold
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

my @parents = ( "F1Mo", "F1Fa", "F0GM" );

my $scaffold_name      = "";
my $in_filename        = "";
my $genotype_filename  = "";
my $tippex_filename    = "";
my $haplotype_filename = "";
my $scaffold_agp_filename = "";

my $log_filename = "";
my $options_okay = GetOptions(
    'scaffold=s'   => \$scaffold_name,
    'input=s'      => \$in_filename,
    'genotypes=s'  => \$genotype_filename,
    'tippex=s'     => \$tippex_filename,
    'haplotypes=s' => \$haplotype_filename,
    'log=s'        => \$log_filename,
    'agp=s'        => \$scaffold_agp_filename,
);

croak
"\nUsage: perl infer_rad_loci_haplotypes.pl -i in_file -g genotypes_file -t tippex_file -a scaffold_agp_filename\n"
  if !$options_okay;

if ( $scaffold_name eq "" ) {
    croak "\nPlease specify an input file with -i\n" if ( $in_filename eq "" );
    croak "\nPlease specify a genotypes file with -g\n"
      if ( $genotype_filename eq "" );
    croak "\nPlease specify a tippex file with -t\n"
      if ( $tippex_filename eq "" );
    croak "\nPlease specify a haploytype output file with -h\n"
      if ( $haplotype_filename eq "" );
    croak "\nPlease specify a log output file with -l\n"
      if ( $log_filename eq "" );
}
else {
    $in_filename =
      "PstI.all.110802.Hmel1-1_primaryScaffolds.bam.$scaffold_name\.tsv";
    $genotype_filename = "genotypes.$scaffold_name\.tsv";
    $tippex_filename   = "tippex.$scaffold_name\.tsv";
    $haplotype_filename =
"PstI.all.110802.Hmel1-1_primaryScaffolds.bam.$scaffold_name\.haplotypes.tsv";
    $log_filename =
"PstI.all.110802.Hmel1-1_primaryScaffolds.bam.$scaffold_name\.haplotypes.txt";
}

croak "\nPlease specify a scaffold AGP file with -a\n" if ($scaffold_agp_filename eq "");

open my $logfile, '>', $log_filename
  or croak "Can't open log output file $log_filename: $OS_ERROR!\n";

open my $hapfile, '>', $haplotype_filename
  or croak "Can't open haplotype output file $haplotype_filename: $OS_ERROR!\n";

#############################################################################
###
### PROCESS AGP FILE TO GET CONTIG BREAKS ON SCAFFOLDS

my %scf_ctg_start;
open my $scaffold_agp, '<', $scaffold_agp_filename
  or croak "Can't open scaffold AGP file $scaffold_agp_filename: $OS_ERROR!\n";
while ( my $scaffold_agp_line = <$scaffold_agp> ) {
    chomp $scaffold_agp_line;
    my @f = split /\t/, $scaffold_agp_line;
    next if ( $f[4] ne "W" );
    $scf_ctg_start{ $f[0] }{ $f[1] } = $f[7];
}
close $scaffold_agp;

#############################################################################
###
### LOAD PREPROCESSED GENOTYPES FOR ALL INDIVIDUALS

my %genotype_range;
my %parent_haplotypes;
open my $genotype_file, '<', $genotype_filename
  or croak "Can't open genotypes file $genotype_filename: $OS_ERROR!\n";
my $genotype_header = <$genotype_file>;
while ( my $genotype_line = <$genotype_file> ) {
    chomp $genotype_line;
    next if ( $genotype_line eq "" );
    my ( $sample, $type, $start, $end, $mother, $father ) = split /\t/,
      $genotype_line;

    if ( $type eq "P" ) {
        $parent_haplotypes{$sample}{$start}{1} = $mother;
        $parent_haplotypes{$sample}{$start}{2} = $father;
    }
    $genotype_range{$sample}{$start} = "$mother$father";
}
close $genotype_file;


#############################################################################
###
### LOAD CLEANING INFORMATION FOR THIS SCAFFOLD

my %tippex_haps;
my %ignore;
my $min;
my $max;

open my $tippex_file, '<', $tippex_filename
  or croak "Can't open tippex file $tippex_filename:$OS_ERROR!\n";

while ( my $tippex_line = <$tippex_file> ) {
    chomp $tippex_line;
    next if ( $tippex_line eq "" );
    my ( $line_type, @f ) = split /\t/, $tippex_line;
    if ( $line_type eq "max" ) {
        $max = $f[0];
    }
    elsif ( $line_type eq "min" ) {
        $min = $f[0];
    }
    elsif ( $line_type eq "ignore" ) {
        $ignore{ $f[0] } = 0;
    }
    elsif ( $line_type eq "hap" ) {
        my ( $loc, $dir, $hap, $allele ) = @f;
        $tippex_haps{$loc}{$dir}{$hap}{$allele} = 0;
    }
}
close $tippex_file;

#############################################################################
###
### LOAD LOCI FROM TSV FILE

my $scf = "";
my %tag;
my %tsv;
my %genotype;

open my $in_file, '<', $in_filename
  or croak "Can't open input file $in_filename: $OS_ERROR!\n";

my $header = <$in_file>;
chomp $header;
print $hapfile
"$header\tAlleleCount\tMoHap\tMoFragLen\tMoSite\tMoNeighbour\tFaHap\tFaFragLen\tFaSite\tFaNeighbour\tSharedFragLen\n";

print $logfile "Loading loci...\n";

while ( my $in_line = <$in_file> ) {
    chomp $in_line;

    my (
        $line_scf,              $loc,
        $dir,                   $allele,
        $seq,                   $mintagfraglen,
        $maxtagfraglen,         $mintaglengc,
        $maxtaglengc,           $totallocusreads,
        $totalallelereads,      $sample,
        $samplelocusreads,      $samplelocusfragments,
        $sampleallelereads,     $sampleallelefragments,
        $sampleallelereadsprop, $sampleallelefragsprop,
        $sampleresfraglen,      $samplematchingtag,
        $minsamplefraglen,      $meansamplefraglen,
        $maxsamplefraglen,      $minsamplelengc,
        $meansamplelengc,       $maxsamplelengc
    ) = split /\t/, $in_line;

    next if ( defined $ignore{$loc} );

    if ( ( $scf ne "" ) && ( $line_scf ne $scf ) ) {
        print $logfile
"More than one scaffold detected: $scf, $line_scf\nThis script only processes one scaffold\n";
        exit;
    }

    $scf = $line_scf;

    if ( ( $loc >= $min ) && ( $loc <= $max ) ) {
        $tsv{$loc}{$dir}{$allele}{$sample} = $in_line;

        map {
            if ( $loc >= $_ )
            {
                $genotype{$loc}{$sample} =
                  ( ( $sample eq "F1Mo" ) || ( $sample eq "F1Fa" ) )
                  ? $genotype_range{$sample}{$_}
                  : get_sample_genotype( $sample, $loc,
                    $genotype_range{$sample}{$_},
                    \%parent_haplotypes );
            }
        } sort { $a <=> $b } keys %{ $genotype_range{$sample} };

        if ( $sampleallelefragments > 0 ) {
            $tag{$loc}{$dir}{$allele}{$sample} = $sampleallelefragments;
        }
    }
}
close $in_file;


#############################################################################
###
### CALCULATE HAPLOTYPES FOR EACH LOCUS AND DELETE INVALID LOCI

my %haplotypes;
my %haplotype_alleles;

print $logfile "Generating haplotypes...\n";
print $logfile "Scaffold\t|| Location\t|| Dir\t||";
print $logfile "\t\tGenotypes\t";
print $logfile "\t||";
print $logfile "\t\tHaplotypes\t";
print $logfile "\t||";
map { print $logfile "\t$_"; } @parents;
print $logfile "\t||";
print $logfile "\n";

foreach my $loc ( sort { $a <=> $b } keys %tag ) {
    foreach my $dir ( keys %{ $tag{$loc} } ) {

        my %genotype_samples;
        my %alt;
        foreach my $sample ( keys %{ $genotype{$loc} } ) {
            if ( ( $sample eq "F1Mo" ) || ( $sample eq "F1Fa" ) ) {
                my @haps = split //, $genotype{$loc}{$sample};

                $haplotypes{$loc}{$dir}{$sample}{mother} = $haps[0];
                $haplotypes{$loc}{$dir}{$sample}{father} = $haps[1];

                $alt{"$sample$haps[0]"} = "$sample$haps[1]";
                $alt{"$sample$haps[1]"} = "$sample$haps[0]";
            }
            elsif ( $sample eq "F0GM" ) {
                next;    # Ignore F0 grandmother
            }
            else {
                $genotype_samples{ $genotype{$loc}{$sample} }++;
            }
        }

        my $haploc_r = $haplotypes{$loc}{$dir};
        my %genotype_allele;
        my $output_string = "";

        $output_string .= "$scf\t|| $loc\t|| $dir\t||";

        # Get genotype for each individual
        foreach my $allele ( keys %{ $tag{$loc}{$dir} } ) {
            foreach my $sample ( keys %{ $tag{$loc}{$dir}{$allele} } ) {
                next if ( $sample =~ /^F/ );    # Skip parents
                next
                  if ( !defined $genotype{$loc}{$sample} )
                  ;    # Skip individuals without genotypes (low coverage)
                $genotype_allele{ $genotype{$loc}{$sample} }{$allele}++;
            }
        }

        # Identify alleles valid across whole genotype groups
        my %valid_alleles;
        my %parental_allele;

        foreach
          my $mohap ( $haploc_r->{"F1Mo"}{mother}, $haploc_r->{"F1Mo"}{father} )
        {
            foreach my $fahap ( $haploc_r->{"F1Fa"}{mother},
                $haploc_r->{"F1Fa"}{father} )
            {
                my $gt = $mohap . $fahap;
                $output_string .= "\t$gt:";

                foreach my $allele (
                    sort { $a <=> $b }
                    keys %{ $genotype_allele{$gt} }
                  )
                {
                    if ( $genotype_allele{$gt}{$allele} ==
                        $genotype_samples{$gt} )
                    {
                        $valid_alleles{$allele}++;
                        $output_string .= "$allele/";

                        map {
                            if (   ( defined( $tag{$loc}{$dir}{$allele}{$_} ) )
                                && ( $tag{$loc}{$dir}{$allele}{$_} > 0 ) )
                            {
                                $parental_allele{$_}{$allele} =
                                  $genotype_allele{$gt}{$allele};
                            }
                        } @parents;

                    }
                    else {
                        delete $genotype_allele{$gt}{$allele};
                    }
                }
                chop $output_string if ( substr( $output_string, -1 ) eq "/" );

            }
        }

        my @genotypes;
        $genotypes[0] =
          $haploc_r->{"F1Mo"}{mother} . $haploc_r->{"F1Fa"}{mother};
        $genotypes[1] =
          $haploc_r->{"F1Mo"}{mother} . $haploc_r->{"F1Fa"}{father};
        $genotypes[2] =
          $haploc_r->{"F1Mo"}{father} . $haploc_r->{"F1Fa"}{mother};
        $genotypes[3] =
          $haploc_r->{"F1Mo"}{father} . $haploc_r->{"F1Fa"}{father};

        my @haplotypes = (
            "F1Mo" . $haploc_r->{"F1Mo"}{mother},
            "F1Mo" . $haploc_r->{"F1Mo"}{father},
            "F1Fa" . $haploc_r->{"F1Fa"}{mother},
            "F1Fa" . $haploc_r->{"F1Fa"}{father}
        );

        # Assign alleles to haplotypes
        foreach my $allele ( keys %valid_alleles ) {
            $haplotype_alleles{ $haplotypes[0] }{$loc}{$dir}{$allele}++
              if ( defined( $genotype_allele{ $genotypes[0] }{$allele} )
                && defined( $genotype_allele{ $genotypes[1] }{$allele} )
                && defined( $parental_allele{"F1Mo"}{$allele} ) );
            $haplotype_alleles{ $haplotypes[1] }{$loc}{$dir}{$allele}++
              if ( defined( $genotype_allele{ $genotypes[2] }{$allele} )
                && defined( $genotype_allele{ $genotypes[3] }{$allele} )
                && defined( $parental_allele{"F1Mo"}{$allele} ) );
            $haplotype_alleles{ $haplotypes[2] }{$loc}{$dir}{$allele}++
              if ( defined( $genotype_allele{ $genotypes[0] }{$allele} )
                && defined( $genotype_allele{ $genotypes[2] }{$allele} )
                && defined( $parental_allele{"F1Fa"}{$allele} ) );
            $haplotype_alleles{ $haplotypes[3] }{$loc}{$dir}{$allele}++
              if ( defined( $genotype_allele{ $genotypes[1] }{$allele} )
                && defined( $genotype_allele{ $genotypes[3] }{$allele} )
                && defined( $parental_allele{"F1Fa"}{$allele} ) );
        }

        $output_string .= "\t||";

        # Fix haps where preset

        foreach my $haplotype (@haplotypes) {
            if ( defined $tippex_haps{$loc}{$dir}{$haplotype} ) {

                # Delete inferred haplotypes
                map { delete $haplotype_alleles{$haplotype}{$loc}{$dir}{$_} }
                  keys %{ $haplotype_alleles{$haplotype}{$loc}{$dir} };

                # Imprint new ones
                map {
                    next if ( $_ eq "-" );
                    $haplotype_alleles{$haplotype}{$loc}{$dir}{$_}++;
                } keys %{ $tippex_haps{$loc}{$dir}{$haplotype} };
            }
        }

        # Clean up haplotypes with less or more than one allele
        my $haplotypes_ok = 0;
        foreach my $haplotype (@haplotypes) {
            $output_string .= "\t$haplotype:";
            if ( keys %{ $haplotype_alleles{$haplotype}{$loc}{$dir} } > 1 ) {
                my $altname = $alt{$haplotype};
                if ( keys %{ $haplotype_alleles{$altname}{$loc}{$dir} } == 1 ) {
                    my $fixed_allele =
                      ( keys %{ $haplotype_alleles{$altname}{$loc}{$dir} } )[0];
                    delete $haplotype_alleles{$haplotype}{$loc}{$dir}
                      {$fixed_allele};
                }
            }
            if ( keys %{ $haplotype_alleles{$haplotype}{$loc}{$dir} } != 1 ) {
                delete $haplotype_alleles{$haplotype}{$loc}{$dir};
            }
            else {
                $haplotypes_ok++;
                foreach my $allele (
                    keys %{ $haplotype_alleles{$haplotype}{$loc}{$dir} } )
                {
                    $output_string .= "$allele/";
                }
                chop $output_string
                  if ( substr( $output_string, -1 ) eq "/" );
            }
        }

        # Check parental genotypes match inferred haplotypes for this locus
        $output_string .= "\t||";
        foreach my $parent (@parents) {
            $output_string .= "\t";
            foreach my $allele ( sort keys %{ $parental_allele{$parent} } ) {
                my @haps =
                  @{ get_sample_haplotypes( $parent, $genotype{$loc}{$parent} )
                  };
                if (
                    (
                        defined(
                            $haplotype_alleles{ $haps[0] }{$loc}{$dir}{$allele}
                        )
                    )
                    || (
                        defined(
                            $haplotype_alleles{ $haps[1] }{$loc}{$dir}{$allele}
                        )
                    )
                  )
                {
                    $output_string .= "$allele/";
                }
            }
            chop $output_string if ( substr( $output_string, -1 ) eq "/" );
        }
        $output_string .= "\t||";

        print $logfile "$output_string";
        print $logfile "\n";
    }
}

#############################################################################
###
### FIND NEIGHBOURING LOCI FOR EACH HAPLOTYPE FOR FRAGMENT LENGTH CALCULATION

my %neighbour;
my %fraglen;
foreach my $sample ( keys %genotype_range ) {

    my %s_haplocs;
    my %s_hapdirs;
    my %s_haps;
    foreach my $loc ( sort { $a <=> $b } keys %genotype ) {
        if ( defined $genotype{$loc}{$sample} ) {
            my ( $mohap, $fahap ) =
              @{ get_sample_haplotypes( $sample, $genotype{$loc}{$sample} ) };
            foreach my $dir ( '+', '-' ) {
                if (   ( defined $haplotype_alleles{$mohap}{$loc} )
                    && ( defined $haplotype_alleles{$mohap}{$loc}{$dir} ) )
                {
                    push @{ $s_haplocs{mother} }, $loc;
                    push @{ $s_hapdirs{mother} }, $dir;
                    push @{ $s_haps{mother} },    $mohap;
                }
                if (   ( defined $haplotype_alleles{$fahap}{$loc} )
                    && ( defined $haplotype_alleles{$fahap}{$loc}{$dir} ) )
                {
                    push @{ $s_haplocs{father} }, $loc;
                    push @{ $s_hapdirs{father} }, $dir;
                    push @{ $s_haps{father} },    $fahap;
                }
            }
        }
    }

    foreach my $hap ( "mother", "father" ) {
        for my $i ( 0 .. $#{ $s_haplocs{$hap} } ) {
            my $dir = $s_hapdirs{$hap}->[$i];
            next
              if (
                keys %{
                    $haplotype_alleles{ $s_haps{$hap}->[$i] }
                      { $s_haplocs{$hap}->[$i] }{$dir}
                } == 0
              );    # Ignore location if no alleles present here

            my $match_dir  = $dir eq '+' ? '-' : '+';
            my $search_inc = $dir eq '+' ? 1   : -1;

            my $candidate_loc = $i + $search_inc;
            my $found_match   = 0;
            while (( $candidate_loc >= 0 )
                && ( $candidate_loc <= $#{ $s_haplocs{$hap} } ) )
            {
                my $chap  = $s_haps{$hap}->[$candidate_loc];
                my $chloc = $s_haplocs{$hap}->[$candidate_loc];
                if (
                    ( !defined $haplotype_alleles{$chap}{$chloc} )
                    || (
                        !defined $haplotype_alleles{$chap}{$chloc}{$match_dir} )
                    || (
                        keys %{ $haplotype_alleles{$chap}{$chloc}{$match_dir} }
                        == 0 )
                    || (   ( $search_inc == 1 )
                        && ( $chloc <= $s_haplocs{$hap}->[$i] + 5 ) )
                    || (   ( $search_inc == -1 )
                        && ( $chloc >= $s_haplocs{$hap}->[$i] - 5 ) )
                  )
                {
                    $candidate_loc += $search_inc;
                    next;
                }
                if ( on_one_contig( $scf, $s_haplocs{$hap}->[$i], $chloc ) ) {
                    $neighbour{$sample}{ $s_haps{$hap}->[$i] }
                      { $s_haplocs{$hap}->[$i] }{$dir} = $chloc;
                    $fraglen{$sample}{ $s_haps{$hap}->[$i] }
                      { $s_haplocs{$hap}->[$i] }{$dir} =
                      ( $dir eq '+' )
                      ? $chloc - $s_haplocs{$hap}->[$i] + 1
                      : $s_haplocs{$hap}->[$i] - $chloc + 1;
                }
                last;
            }
        }
    }
}

#############################################################################
###
### OUTPUT LOCI WITH VALID HAPLOTYPES TO TSV FILE WITH HAPLOTYPE INFORMATION

foreach my $loc ( sort { $a <=> $b } keys %tag ) {

    foreach my $dir ( sort keys %{ $tag{$loc} } ) {

        foreach my $allele ( sort { $a <=> $b } keys %{ $tag{$loc}{$dir} } ) {

            my %fraglen_allele;
            foreach my $sample ( keys %{ $genotype{$loc} } ) {
                foreach my $haplotype (
                    @{
                        get_sample_haplotypes( $sample,
                            $genotype{$loc}{$sample} )
                    }
                  )
                {
                    if (
                        (
                            defined $haplotype_alleles{$haplotype}{$loc}{$dir}
                            {$allele}
                        )
                        && ( defined $fraglen{$sample}{$haplotype}{$loc}{$dir} )
                      )
                    {
                        $fraglen_allele{ $fraglen{$sample}{$haplotype}{$loc}
                              {$dir} }++;
                    }
                }
            }

            # Create sample set for this locus,
            # sorted by genotype, with parents at the start
            my @samples;
            map {
                if ( $_ !~ /^F/ ) { push @samples, $_; }
            } keys %{ $genotype{$loc} };
            @samples =
              sort { $genotype{$loc}{$a} cmp $genotype{$loc}{$b} || $a <=> $b }
              @samples;
            unshift @samples, ( "F0GM", "F1Mo", "F1Fa" );

            foreach my $sample (@samples) {
                my @sample_haplotypes =
                  @{ get_sample_haplotypes( $sample, $genotype{$loc}{$sample} )
                  };

                my $allele_count = 0;
                foreach my $haplotype (@sample_haplotypes) {
                    $allele_count++
                      if (
                        defined(
                            $haplotype_alleles{$haplotype}{$loc}{$dir}{$allele}
                        )
                      );
                }
                print $hapfile
                  "$tsv{$loc}{$dir}{$allele}{$sample}\t$allele_count";
                my %fraglens;
                if (
                    (
                        defined $haplotype_alleles{ $sample_haplotypes[0] }
                        {$loc}{$dir}{$allele}
                    )
                    && (
                        defined $fraglen{$sample}{ $sample_haplotypes[0] }{$loc}
                        {$dir} )

 #                    && (
 #                        on_one_contig(
 #                            $scf,
 #                            $loc,
 #                            $neighbour{$sample}{ $sample_haplotypes[0] }{$loc}
 #                              {$dir}
 #                        )
 #                    )
                  )
                {
                    print $hapfile
"\t$sample_haplotypes[0]\t$fraglen{$sample}{$sample_haplotypes[0]}{$loc}{$dir}\t$loc\t$neighbour{$sample}{$sample_haplotypes[0]}{$loc}{$dir}";

                    $fraglens{ $fraglen{$sample}{ $sample_haplotypes[0] }{$loc}
                          {$dir} }++;
                }
                else {
                    print $hapfile "\tNA\tNA\tNA\tNA";
                }
                if (
                    (
                        defined $haplotype_alleles{ $sample_haplotypes[1] }
                        {$loc}{$dir}{$allele}
                    )
                    && (
                        defined $fraglen{$sample}{ $sample_haplotypes[1] }{$loc}
                        {$dir} )

 #                    && (
 #                        on_one_contig(
 #                            $scf,
 #                            $loc,
 #                            $neighbour{$sample}{ $sample_haplotypes[1] }{$loc}
 #                              {$dir}
 #                        )
 #                    )
                  )
                {
                    print $hapfile
"\t$sample_haplotypes[1]\t$fraglen{$sample}{$sample_haplotypes[1]}{$loc}{$dir}\t$loc\t$neighbour{$sample}{$sample_haplotypes[1]}{$loc}{$dir}";

                    $fraglens{ $fraglen{$sample}{ $sample_haplotypes[1] }{$loc}
                          {$dir} }++;
                }
                else {
                    print $hapfile "\tNA\tNA\tNA\tNA";
                }
                print $hapfile "\t";

                if ( keys %fraglens == 1 ) {
                    my $fraglen = ( keys %fraglens )[0];
                    print $hapfile $fraglen;
                }
                elsif ( keys %fraglen_allele == 1 ) {
                    my $fraglen = ( keys %fraglen_allele )[0];
                    print $hapfile $fraglen;
                }
                else {
                    print $hapfile "NA";

                }
                print $hapfile "\n";
            }
        }
    }
}

close $logfile;
close $hapfile;



#############################################################################
###
### SUBROUTINES

sub get_sample_genotype {
    my ( $sample, $loc, $sample_haps, $parent_haps_ref ) = @_;

    my @sample_haps = split //, $sample_haps;

    my $genotype;
    my %ploc;
    foreach my $parent ( "F1Mo", "F1Fa" ) {
        map {
            if ( $loc >= $_ )
            {
                $ploc{$parent} = $parent_haps_ref->{$parent}{$_};
            }
        } sort { $a <=> $b } keys %{ $parent_haps_ref->{$parent} };
    }
    $genotype =
      $ploc{"F1Mo"}->{ $sample_haps[0] } . $ploc{"F1Fa"}->{ $sample_haps[1] };
    return $genotype;
}

sub get_sample_haplotypes {
    my ( $sample, $genotypes ) = @_;
    my @gt = split //, $genotypes;
    my @haplotypes;
    if ( $sample eq "F1Mo" ) {
        push @haplotypes, "F1Mo$gt[0]";
        push @haplotypes, "F1Mo$gt[1]";

    }
    elsif ( $sample eq "F1Fa" ) {
        push @haplotypes, "F1Fa$gt[0]";
        push @haplotypes, "F1Fa$gt[1]";

    }
    else {
        push @haplotypes, "F1Mo$gt[0]";
        push @haplotypes, "F1Fa$gt[1]";
    }

    return \@haplotypes;
}

sub on_one_contig {
    my ( $scf, $a, $b ) = @_;
    my $min = min( $a, $b );
    my $max = max( $a, $b );
    foreach my $start ( keys %{ $scf_ctg_start{$scf} } ) {
        return 0 if ( ( $min <= $start ) && ( $start <= $max ) );
    }
    return 1;
}
