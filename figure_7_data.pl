#!/usr/bin/env perl

# figure_7_data.pl
# John W Davey johnomics@gmail.com

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

my @genotypes = ( "AB", "AC", "BC", "CC" );
my %parents = (
    "F1Mo" => "AC",
    "F1Fa" => "BC",
    "F0GM" => "AB",
);

my $tsv_filename = "";

my $stacks_dir  = "";
my $radtags_dir = "";

my $vcf_filename = "";

my $radmarkers_filename = "";
my $catalog_filename    = "";
my $haplotypes_filename = "";

my $options_okay = GetOptions(
    'tsv=s'        => \$tsv_filename,
    'stacks=s'     => \$stacks_dir,
    'radtags=s'    => \$radtags_dir,
    'vcf=s'        => \$vcf_filename,
    'markers=s'    => \$radmarkers_filename,
    'haplotypes=s' => \$haplotypes_filename,
    'catalog=s'    => \$catalog_filename,
);

croak
"\nUsage: perl rad_compare_tool.pl -t tsv_file -s stacks_dir -r radtags_dir -v vcf_filename -h stacks_haplotypes_file -m radmarkers_file -c catalog_file\n"
  if !$options_okay;

croak "Please specify a TSV file with -t\n" if ( $tsv_filename eq "" );

croak "Please specify a directory of Stacks output with -s\n"
  if ( $stacks_dir eq "" );
croak "Please specify a directory of RADtags output with -r\n"
  if ( $radtags_dir eq "" );

croak "Please specify a GATK VCF file with -v\n" if ( $vcf_filename eq "" );
croak "Please specify a RADmarkers output file with -m"
  if ( $radmarkers_filename eq "" );
croak "Please specify a Stacks haplotypes file with -h"
  if ( $haplotypes_filename eq "" );
croak "Please specify a Stacks catalog file with -c"
  if ( $catalog_filename eq "" );

#############################################################################
###
### LOAD HETEROZYGOUS CALLS FROM TSV FILE

my %rad_gts;
my %rad_haps;
my %sample;
my %seq_lookup;

print STDERR "Loading TSV file...\n";

open my $tsv_file, '<', $tsv_filename
  or croak "Can't open $tsv_filename: $OS_ERROR!\n";
my $tsv_header = <$tsv_file>;
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
    next if $f[3] > 3;    #Ignore alleles 4 and above
    my %f;
    for my $i ( 0 .. $#f ) {
        $f{ $tsv_field{$i} } = $f[$i];
    }

    next if ( $f{Sample} eq "40" );

    # Get parental genotype
    if ( $f{Sample} =~ /^F/ ) {
        $rad_gts{ $f{Chr} }{ $f{Loc} }{ $f{Dir} }{ $f{Sample} }{gt} =
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
        $rad_gts{ $f{Chr} }{ $f{Loc} }{ $f{Dir} }{ $f{Sample} }{gt} =
          $genotypes[$genotype_index];
    }

    # Initialise haplotypes to null based on genotype
    if ( $f{Allele} == $first_allele ) {
        my @haps = split //,
          $rad_gts{ $f{Chr} }{ $f{Loc} }{ $f{Dir} }{ $f{Sample} }{gt};
        map {
            $rad_haps{ $f{Chr} }{ $f{Loc} }{ $f{Dir} }{$_}{allele}  = "-";
            $rad_haps{ $f{Chr} }{ $f{Loc} }{ $f{Dir} }{$_}{fraglen} = 0;
        } @haps;
    }

    # Add sample name to hash of samples for later use
    $sample{ $f{Sample} }++;

    if ( !defined $seq_lookup{ $f{Seq} }
        && ( $f{MoHap} ne "NA" || $f{FaHap} ne "NA" ) )
    {
        $seq_lookup{ $f{Seq} }{scf}    = $f{'Chr'};
        $seq_lookup{ $f{Seq} }{loc}    = $f{'Loc'};
        $seq_lookup{ $f{Seq} }{dir}    = $f{'Dir'};
        $seq_lookup{ $f{Seq} }{allele} = $f{'Allele'};
    }

    # Load real haplotypes where present
    for my $parhap ( $f{MoHap}, $f{FaHap} ) {
        next if ( $parhap eq "NA" );
        my $hap = substr $parhap, -1, 1;
        $rad_haps{ $f{Chr} }{ $f{Loc} }{ $f{Dir} }{$hap}{allele} = $f{Allele};
        $rad_haps{ $f{Chr} }{ $f{Loc} }{ $f{Dir} }{$hap}{fraglen} =
          $parhap eq $f{MoHap} ? $f{MoFragLen} : $f{FaFragLen};

        next if ($f{AlleleCount} != 1);
        $rad_haps{ $f{Chr} }{ $f{Loc} }{ $f{Dir} }{$hap}{reads}{ $f{Sample} } =
          $f{SampleAlleleReads};
    }
}
close $tsv_file;

#############################################################################
###
### DISCARD GENOTYPES WITH NULL HAPLOTYPES AND MISSING ALLELES

my %rad_hets;
my %rad_position;
my %variants;
foreach my $scf ( keys %rad_gts ) {
    foreach my $loc ( keys %{ $rad_gts{$scf} } ) {
        foreach my $dir ( keys %{ $rad_gts{$scf}{$loc} } ) {

            map {
                next if $rad_haps{$scf}{$loc}{$dir}{$_}{allele} eq "-";
            } (keys %{$rad_haps{$scf}{$loc}{$dir}});

            # Make list of positions for this locus based on locus direction
            my $dir_mult = $dir eq '+' ? 1 : -1;
            my $read1end = $loc + 150 * $dir_mult;
            my @positions =
                $loc < $read1end
              ? $loc + 1 .. $read1end
              : $read1end .. $loc - 1;
            map { $rad_position{$scf}{$_}{$loc} = 0 } @positions;

            $variants{$scf}{$loc}{$dir}{snps} = 0;
            $variants{$scf}{$loc}{$dir}{indels} = 0;

            foreach my $sample ( keys %{ $rad_gts{$scf}{$loc}{$dir} } ) {
                my @haps = split //, $rad_gts{$scf}{$loc}{$dir}{$sample}{gt};
                next
                  if missing_alleles( $rad_haps{$scf}{$loc}{$dir},
                    \@haps, $sample );

                my $a1 = $rad_haps{$scf}{$loc}{$dir}{ $haps[0] }{allele};
                my $a2 = $rad_haps{$scf}{$loc}{$dir}{ $haps[1] }{allele};
                if ( $a1 ne "-" && $a2 ne "-" && $a1 ne $a2 ) {
                    $rad_hets{$scf}{$loc}{$dir}{$sample}++;
                }
            }
        }
    }
}

sub missing_alleles {
    my ( $rad_hap_ref, $haps_ref, $sample ) = @_;
    return (
        (
            defined $rad_hap_ref->{ $haps_ref->[0] }{reads}{$sample}
              && $rad_hap_ref->{ $haps_ref->[0] }{reads}{$sample} == 0
        )
          || ( defined $rad_hap_ref->{ $haps_ref->[1] }{reads}{$sample}
            && $rad_hap_ref->{ $haps_ref->[1] }{reads}{$sample} == 0 )
    );
}

#############################################################################
###
### LOAD STACKS OUTPUT FOR ALL INDIVIDUALS

print STDERR "Loading Stacks output...\n";

my %stacks_loci;
foreach my $sample ( keys %sample ) {
    my $stacks_filestub = "$stacks_dir/$sample\_1";

    my $stacks_tags_filename = "$stacks_filestub\.tags.tsv";
    open my $stacks_tags_file, '<', $stacks_tags_filename
      or croak "Can't open Stacks tags file $stacks_tags_filename:$OS_ERROR\n";
    while ( my $stacks_tag_line = <$stacks_tags_file> ) {
        chomp $stacks_tag_line;
        if ( $stacks_tag_line =~ /primary/ ) {
            my @f = split /\t/, $stacks_tag_line;
            my $seq = $f[9];
            $stacks_loci{$sample}{ $f[7] }{ $f[2] } = $seq;
            if ( defined $seq_lookup{$seq} ) {
                $rad_gts{ $seq_lookup{$seq}{scf} }{ $seq_lookup{$seq}{loc} }
                  { $seq_lookup{$seq}{dir} }{$sample}{stacks}{ $f[2] }
                  { $f[7] } = $seq_lookup{$seq}{allele};
            }
        }
    }
    close $stacks_tags_file;
}

#############################################################################
###
### LOAD STACKS CATALOGUE

my %catalog;
open my $catalog_file, '<', $catalog_filename
  or croak "Can't open Stacks catalog file $catalog_filename:$OS_ERROR\n";
while ( my $catalog_line = <$catalog_file> ) {
    chomp $catalog_line;
    my @f = split /\t/, $catalog_line;
    my $catalog_id = $f[2];
    foreach my $allele ( split /,/, $f[8] ) {
        my ( $parent_id, $allele_id ) = split /_/, $allele;
        my $parent = $parent_id == 1 ? "F1Fa" : $parent_id == 2 ? "F1Mo" : "-";
        $catalog{$parent}{$allele_id} = $catalog_id;
    }
}
close $catalog_file;

#############################################################################
###
### LOAD STACKS GENOTYPES AND LINK TO SAMPLES

my %stacks_haplotypes;
open my $haplotypes_file, '<', $haplotypes_filename
  or croak "Can't open Stacks haplotypes file $haplotypes_filename:$OS_ERROR\n";
my $haplotypes_header = <$haplotypes_file>;
chomp $haplotypes_header;
my @haplotypes_fields = split /\t/, $haplotypes_header;
my %haplotypes_field;
for my $i ( 0 .. $#haplotypes_fields ) {
    $haplotypes_field{$i} = $haplotypes_fields[$i];
}

while ( my $haplotypes_line = <$haplotypes_file> ) {
    chomp $haplotypes_line;
    my @f = split /\t/, $haplotypes_line;
    my %f;
    for my $i ( 0 .. $#f ) {
        $f{ $haplotypes_field{$i} } = $f[$i];
    }

    foreach my $sample ( keys %f ) {
        next if ( $sample !~ /_1$/ );
        $sample =~ s/_1//;
        $stacks_haplotypes{$sample}{ $f{"Catalog ID"} } = $f{"$sample\_1"};
    }
}
close $haplotypes_file;

#############################################################################
###
### LOAD RADMARKERS OUTPUT

print STDERR "Loading RADmarkers output...\n";

open my $radmarkers_file, '<', $radmarkers_filename
  or croak "Can't open RADmarkers file $radmarkers_filename:$OS_ERROR\n";
my $radmarkers_header = <$radmarkers_file>;
chomp $radmarkers_header;
my @radmarkers_fields = split /\t/, $radmarkers_header;
my %radmarkers_field;
for my $i ( 0 .. $#radmarkers_fields ) {
    $radmarkers_field{$i} = $radmarkers_fields[$i];
}

while ( my $radallele = <$radmarkers_file> ) {
    chomp $radallele;
    my @f = split /\t/, $radallele;
    my %f;
    for my $i ( 0 .. $#f ) {
        $f{ $radmarkers_field{$i} } = $f[$i];
    }
    if ( defined $seq_lookup{ $f{Tag} } ) {
        foreach my $sample ( keys %sample ) {
            push @{ $rad_gts{ $seq_lookup{ $f{Tag} }{scf} }
                  { $seq_lookup{ $f{Tag} }{loc} }{ $seq_lookup{ $f{Tag} }{dir} }
                  {$sample}{radtools}{ $f{ClusterID} } }, $f{$sample};
        }
    }
}
close $radmarkers_file;



#############################################################################
###
### LOAD GATK VCF FILE

my %homhet = (
    "./." => "./.",
    "0/0" => "hom",
    "0/1" => "het",
    "0/2" => "het",
    "0/3" => "het",
    "1/1" => "hom",
    "1/2" => "het",
    "2/2" => "hom",
    "1/3" => "het",
    "2/3" => "het",
    "3/3" => "hom",
);

print STDERR "Loading GATK output...\n";
open my $vcf_file, '<', $vcf_filename
  or croak "Can't open $vcf_filename: $OS_ERROR!\n";

my %vcf_field;

VCF_LINE:
while ( my $vcf_line = <$vcf_file> ) {
    chomp $vcf_line;

    # Parse header line to get field names (including sample names)
    if ( $vcf_line =~ /#CHROM/ ) {
        $vcf_line = substr $vcf_line, 1;    # Strip '#' from CHROM field
        my @f = split /\t/, $vcf_line;
        for my $i ( 0 .. $#f ) {
            if ( ( $i > 8 ) && ( $f[$i] =~ /PstI\./ ) ) {
                $f[$i] = substr $f[$i], 5;
            }
            $vcf_field{$i} = $f[$i];
        }
        next;
    }

    # Skip header comments
    next if ( $vcf_line =~ /^##/ );

    # Process scaffold position on VCF line
    my @f = split /\t/, $vcf_line;
    my %f;
    my %format_field;
    my $locus;
    my $dir;
    my %tsv_vcf;
    my $variant = 0;
    for my $i ( 0 .. $#f ) {
        $f{ $vcf_field{$i} } = $f[$i];

        # Skip paired end bases
        if ( $vcf_field{$i} eq "POS" ) {
            next VCF_LINE
              if ( !defined $rad_position{ $f{CHROM} }{ $f{POS} } );
        }

        # Check if field is variant
        $variant++ if ( ( $vcf_field{$i} eq "ALT" ) && ( $f[$i] ne '.' ) );
        $variant++
          if ( ( $vcf_field{$i} eq "INFO" ) && ( $f[$i] =~ /^INDEL/ ) );

        # Load format field positions (GQ, PL etc)
        if ( $vcf_field{$i} eq "FORMAT" ) {
            my @format_fields = split /:/, $f[$i];
            for my $i ( 0 .. $#format_fields ) {
                $format_field{ $format_fields[$i] } = $i;
            }
        }

        # Count SNPs and indels
        if ( $vcf_field{$i} eq "ALT" ) {
            foreach
              my $locus ( keys %{ $rad_position{ $f{CHROM} }{ $f{POS} } } )
            {
                my $dir = $locus < $f{POS} ? "+" : "-";
                if ( $f{ALT} ne '.' ) {
                    my $primary_alt = ( split /,/, $f{ALT} )[0];
                    if ( length( $f{REF} ) == length($primary_alt) ) {
                        $variants{ $f{CHROM} }{$locus}{$dir}{snps}++;
                    }
                    else {
                        $variants{ $f{CHROM} }{$locus}{$dir}{indels}++;
                    }
                }
            }
        }

        # Load sample genotypes
        if ( ( $i > 8 ) && ( defined $sample{ $vcf_field{$i} } ) ) {
            my @sample_fields = split /:/, $f[$i];
            my $vcf_genotype = $sample_fields[ $format_field{GT} ];

            foreach
              my $locus ( keys %{ $rad_position{ $f{CHROM} }{ $f{POS} } } )
            {
                my $dir = $locus < $f{POS} ? "+" : "-";
                $rad_gts{ $f{CHROM} }{$locus}{$dir}{ $vcf_field{$i} }{gatk_gt}{$homhet{$vcf_genotype}}++;
            }
        }
    }

}
close $vcf_file;


#############################################################################
###
### VALIDATE HETEROZYGOUS GENOTYPES FROM EACH TOOL AND OUTPUT LOCI INFO

print STDERR "Validating heterozygous genotypes...\n";

my %fraglen_gt;

foreach my $scf ( keys %rad_hets ) {
    foreach my $loc ( keys %{ $rad_hets{$scf} } ) {
        foreach my $dir ( keys %{ $rad_hets{$scf}{$loc} } ) {

            # Get Stacks catalogue ID for this locus from parents
            my %locus_catalog_ids;
            foreach my $parent ( "F1Mo", "F1Fa" ) {
                next
                  if ( !defined $rad_gts{$scf}{$loc}{$dir}{$parent}{stacks} );
                foreach my $stacks_locus (
                    keys %{ $rad_gts{$scf}{$loc}{$dir}{$parent}{stacks} } )
                {
                    if ( defined $catalog{$parent}{$stacks_locus} ) {
                        $locus_catalog_ids{ $catalog{$parent}
                              {$stacks_locus} }++;
                    }
                }
            }

            foreach my $sample ( keys %{ $rad_hets{$scf}{$loc}{$dir} } ) {

# Assign Stacks genotypes for each sample genotype based on matching catalog ID(s)
                $rad_gts{$scf}{$loc}{$dir}{$sample}{stacks_gt} = "--";
                if ( keys %locus_catalog_ids > 1 ) {
                    $rad_gts{$scf}{$loc}{$dir}{$sample}{stacks_gt} = "NC";
                }
                elsif ( keys %locus_catalog_ids == 0 ) {
                    $rad_gts{$scf}{$loc}{$dir}{$sample}{stacks_gt} = "NA";
                }
                else {
                    my $catalog_id = ( keys %locus_catalog_ids )[0];
                    $rad_gts{$scf}{$loc}{$dir}{$sample}{stacks_gt} =
                      $stacks_haplotypes{$sample}{$catalog_id};
                }

                $rad_gts{$scf}{$loc}{$dir}{$sample}{radtools_gt} = "--";
                if ( !defined $rad_gts{$scf}{$loc}{$dir}{$sample}{radtools} ) {
                    $rad_gts{$scf}{$loc}{$dir}{$sample}{radtools_gt} = "NA";
                }
                elsif (
                    keys %{ $rad_gts{$scf}{$loc}{$dir}{$sample}{radtools} } >
                    1 )
                {
                    $rad_gts{$scf}{$loc}{$dir}{$sample}{radtools_gt} = "NC";
                }
                else {
                    my $rt_clusterid = (
                        keys %{ $rad_gts{$scf}{$loc}{$dir}{$sample}{radtools} }
                    )[0];
                    my @rt_alleles;
                    map {
                        if ( $_ ne "NA" ) { push @rt_alleles, $_ }
                      } @{ $rad_gts{$scf}{$loc}{$dir}{$sample}{radtools}
                          {$rt_clusterid} };
                    $rad_gts{$scf}{$loc}{$dir}{$sample}{radtools_gt} = join "/",
                      @rt_alleles;
                }
                my @haps = split //, $rad_gts{$scf}{$loc}{$dir}{$sample}{gt};
                my $fraglena = $rad_haps{$scf}{$loc}{$dir}{$haps[0]}{fraglen};
                my $fraglenb = $rad_haps{$scf}{$loc}{$dir}{$haps[1]}{fraglen};
                my $hapa = $haps[0];
                my $hapb = $haps[1];
                if ($fraglenb > $fraglena) {
                    my $temp_fraglen = $fraglenb;
                    $fraglenb = $fraglena;
                    $fraglena = $temp_fraglen;
                    $hapa = $haps[1];
                    $hapb = $haps[0];
                }
                
               $rad_gts{$scf}{$loc}{$dir}{$sample}{readprop} = sprintf "%3.2f",  $rad_haps{$scf}{$loc}{$dir}{$hapa}{reads}{$sample}/($rad_haps{$scf}{$loc}{$dir}{$hapa}{reads}{$sample}+$rad_haps{$scf}{$loc}{$dir}{$hapb}{reads}{$sample}); 

                print STDERR "$scf\t$loc\t$dir\t$sample\t$fraglena\t$fraglenb\t$variants{$scf}{$loc}{$dir}{snps}\t$variants{$scf}{$loc}{$dir}{indels}\t$rad_gts{$scf}{$loc}{$dir}{$sample}{readprop}\t$rad_gts{$scf}{$loc}{$dir}{$sample}{stacks_gt}\t$rad_gts{$scf}{$loc}{$dir}{$sample}{radtools_gt}";
                map {
                    print STDERR "\t$_:$rad_gts{$scf}{$loc}{$dir}{$sample}{gatk_gt}{$_}"
                } sort keys %{$rad_gts{$scf}{$loc}{$dir}{$sample}{gatk_gt}};
                print STDERR "\n";

                next if ($variants{$scf}{$loc}{$dir}{indels} > 0);
                next if ($rad_gts{$scf}{$loc}{$dir}{$sample}{stacks_gt} eq "NC" && $rad_gts{$scf}{$loc}{$dir}{$sample}{radtools_gt} eq "NC");

                if (!defined $fraglen_gt{$fraglena}{$fraglenb}) {
                    $fraglen_gt{$fraglena}{$fraglenb}{total} = 0;
                    $fraglen_gt{$fraglena}{$fraglenb}{Stacks} = 0;
                    $fraglen_gt{$fraglena}{$fraglenb}{RADtools} = 0;
                    $fraglen_gt{$fraglena}{$fraglenb}{GATK} = 0;
                }
                
                $fraglen_gt{$fraglena}{$fraglenb}{total}++;

                push @{$fraglen_gt{$fraglena}{$fraglenb}{readprop}}, $rad_gts{$scf}{$loc}{$dir}{$sample}{readprop};

                $fraglen_gt{$fraglena}{$fraglenb}{"Stacks Clusters"}++ if ($rad_gts{$scf}{$loc}{$dir}{$sample}{stacks_gt} ne "NC");

                my @stacks_alleles = split "/", $rad_gts{$scf}{$loc}{$dir}{$sample}{stacks_gt};
                $fraglen_gt{$fraglena}{$fraglenb}{"Stacks Genotypes"}++ if (@stacks_alleles == 2);
                
                my @radtools_alleles = split "/", $rad_gts{$scf}{$loc}{$dir}{$sample}{radtools_gt};
                $fraglen_gt{$fraglena}{$fraglenb}{"RADtools Clusters"}++ if (@radtools_alleles == 2);

                $fraglen_gt{$fraglena}{$fraglenb}{GATK}++ if (defined $rad_gts{$scf}{$loc}{$dir}{$sample}{gatk_gt}{het});
            }
        }
    }
}


#############################################################################
###
### OUTPUT TOOL PERFORMANCE ON EACH FRAGMENT LENGTH PAIR

print "FragLen1\tFragLen2\tTool\tProportion\n";
foreach my $fraglena (sort {$a<=>$b} keys %fraglen_gt) {
    foreach my $fraglenb (sort {$a<=>$b} keys %{$fraglen_gt{$fraglena}}) {
        my $av_readprop = sprintf "%3.2f", sum(@{$fraglen_gt{$fraglena}{$fraglenb}{readprop}})/@{$fraglen_gt{$fraglena}{$fraglenb}{readprop}};
        print "$fraglena\t$fraglenb\tReads\t$av_readprop\n";
        foreach my $tool ("Stacks Clusters", "Stacks Genotypes", "RADtools Clusters", "GATK") {
            $fraglen_gt{$fraglena}{$fraglenb}{$tool} = 0 if (!defined $fraglen_gt{$fraglena}{$fraglenb}{$tool});
            my $tool_ratio = sprintf "%3.2f", 1 - $fraglen_gt{$fraglena}{$fraglenb}{$tool} / $fraglen_gt{$fraglena}{$fraglenb}{total};
            print "$fraglena\t$fraglenb\t$tool\t$tool_ratio\n";

        }
    }
}
