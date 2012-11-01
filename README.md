# Special features of RAD Sequencing data

This document describes a set of scripts and commands to replicate the figures in the following paper, and how to use the scripts on other data sets. Please refer to the manuscript for further explanation and justification for the methods used.

Davey JW, Cezard T, Fuentes-Utrilla P, Eland C, Gharbi K, Blaxter M. Special features of RAD Sequencing data: implications for genotyping. [Molecular Ecology, Early View](http://onlinelibrary.wiley.com/doi/10.1111/mec.12084/full).

This document is primarily provided for transparency and many scripts are specific to the analyses reported in the paper. However, a good faith attempt has been made to generalize several tools that may be of use on other data sets, notably `simulate_RAD_fragments.pl`, `rad_locus_perfect_matches.pl`, `summarise_rad_loci.pl`, and `infer_rad_loci_per_scaffold_from_bam.pl`. Be warned that these scripts have not been tested on any other data sets or machines, so errors are likely to occur. If you believe you have found a bug in these scripts, please contact [John Davey](mailto:johnomics@gmail.com) or make a pull request to the [github repository](https://github.com/johnomics/special_features_of_rad_sequencing).

These scripts are uploaded to Dryad repository DOI [10.5061/dryad.218p2](http://datadryad.org/resource/doi:10.5061/dryad.218p2).

## FILE LIST 

The following files are included in this distribution. These files are used for the commands below. Novel file formats are described at the end of this document. Intermediate files created by the commands listed below are not included, except for the files describing the C. elegans and Heliconius RAD loci.

### Scripts


````
figure_6_data.pl - generate missing data statistics for Figure 6 of the paper
figure_7_data.pl - generate tool comparison statistics for Figure 7 of the paper
infer_rad_loci_haplotypes.pl - calculate haplotypes for a Heliconius PstI cross based on inferred RAD loci
infer_rad_loci_per_scaffold_from_bam.pl - infer the positions of RAD loci from a full alignment of raw reads to a reference genome
infer_rad_loci_statistics.pl - calculate summary statistics for inferred RAD loci based on aligned reads
plot_restriction_sites_fraglen.R - plot RAD loci and other features along a scaffold with aligned reads
simulate_rad_fragments.pl - generate a set of in silico restriction fragment sequences, given a reference genome and a restriction enzyme cut site
summarise_rad_loci.pl - calculate summary statistics for the in silico RAD loci based on aligned reads
figure_6_plot.R - plot Figure 6 of the paper
figure_7_plot.R - plot Figure 7 of the paper
rad_locus_perfect_matches.pl - generate a summary of sequence reads aligned to in silico RAD loci
````

### Data files

````
celegans_rad_samples.pools - sample names and barcodes for the C. elegans PstI RAD library
toolcomp.fastqfiles.pools - sample names and barcodes for the Heliconius PstI RAD library, for running RADtools on Heliconius data

chr18.scaffolds.clean.names.tsv - list of Heliconius Chromosome 18 scaffolds with usable RAD loci

psti.rad.read.counts.tsv - read counts per sample for Heliconius PstI RAD library

chr18.scaffolds.clean.loci.tsv - cleaned inferred Heliconius RAD loci per sample with many summary statistics
c_elegans.loci.tsv - list of in silico C. elegans RAD loci with summary statistics for perfectly aligned PstI RAD reads
c_elegans.WS229.genomic.PstI.fragments.txt - in silico C. elegans PstI fragments, for reference to c_elegans.loci.tsv

heliconius_scaffold_files - folder containing genotypes, tippex files, and PDF plots for all Heliconius chromosome 18 scaffolds with usable RAD loci
````

Additional data files including the raw C. elegans data `111206_0150_BD0BP3ACXX_1_1.sanfastq.gz` and `111206_0150_BD0BP3ACXX_1_2.sanfastq.gz` and the Heliconius BAM file `PstI.all.110802.Hmel1-1_primaryScaffolds.chr18.bam` and the MiSeq sequence lanes for Figure S3 are available under European Nucleotide Archive project ERP001757.


## C. ELEGANS DATA SET

### Figure 1

Raw C. elegans reads were demultiplexed with [RADpools](http://radseq.info) using the barcode file `celegans_rad_samples.pools`:

```
RADpools -t 95 -i <(gunzip -c 111206_0150_BD0BP3ACXX_1_1.sanfastq.gz)  -p <(gunzip -c 111206_0150_BD0BP3ACXX_1_2.sanfastq.gz) -m 8 -f -d celegans_rad_samples -e TGCAG -o -s
```

Demultiplexed FASTQ files are available in the European Nucleotide Archive, project ERP001757.

Raw reads for N2_14_1 sample aligned to C. elegans genome version [WS229](http://wiki.wormbase.org/index.php/WS229) using BWA [v0.6.1](http://sourceforge.net/projects/bio-bwa/files/bwa-0.6.1.tar.bz2/download), SAMtools [v0.1.17](http://sourceforge.net/projects/samtools/files/samtools/0.1.17/) and Picard tools [v1.55](http://sourceforge.net/projects/picard/files/picard-tools/1.55/).

```
bwa aln  -t 4 c_elegans.WS229.genomic.fa N2_14_1_1.fastq > N2_14_1_1.sai
bwa aln  -t 4 c_elegans.WS229.genomic.fa N2_14_1_2.fastq > N2_14_1_2.sai
bwa sampe -r "@RG\tID:N2_14_1\tLB:N2_14_1\tPL:ILLUMINA\tSM:N2_14_1" c_elegans.WS229.genomic.fa N2_14_1_1.sai N2_14_1_2.sai N2_14_1_1.fastq N2_14_1_2.fastq | samtools view -bS - > N2_14_1.bam
java -Xmx4G -jar SortSam.jar I=N2_14_1.bam O=N2_14_1_sorted.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT
```

BAM and genome viewed using IGV [v2.1.16](http://www.broadinstitute.org/igv/node/236), region `CHROMOSOME_I:6,649,884-6,658,814`.

### Figure 2

A quick Unix hack was used to select reads perfectly matching N2_14_1 barcode and generate histogram from read counts:

```
gunzip -c 111206_0150_BD0BP3ACXX_1_1.sanfastq.gz | perl -ne 'print if (/^ACCATTGCAG([ACGTN]+)$/)' | cut -c11- | sort | uniq -c | awk '{print $1}' | sort -n | uniq -c > 111206_0150_BD0BP3ACXX_1_1.N2_14_1.read.counts.txt
```

```r
> library(ggplot2)
> library(scales)
> read.csv("111206_0150_BD0BP3ACXX_1_1.N2_14_1.read.counts.txt",header=F,sep="")->N2.14.1.read.hist
> colnames(N2.14.1.read.hist)<-c("Uniques","Reads")
> ggplot(N2.14.1.read.hist,aes(Reads,Uniques))+geom_point()+scale_x_continuous(name="Read Depth", limits=c(0,1000))+ scale_y_log10(labels=comma_format(),name="Unique Sequences")+theme_bw(base_family="Palatino")
> dev.off()
```


### Figure 3

Generate PstI RAD fragments from C. elegans genome:
```
simulate_rad_fragments.pl -i c_elegans.WS229.genomic.fa -e CTGCAG > c_elegans.WS229.genomic.PstI.fragments.txt
```

Calculate coverage of each RAD locus by perfectly matching reads:

```
rad_locus_perfect_matches.pl -f c_elegans.WS229.genomic.PstI.fragments.txt -l 96 -e C*TGCA*G -1 <(gunzip -c 111206_0150_BD0BP3ACXX_1_1.sanfastq.gz) -2 <(gunzip -c 111206_0150_BD0BP3ACXX_1_2.sanfastq.gz) -b celegans_rad_samples.pools -c 111206_0150_BD0BP3ACXX_1.locus.pileup.csv -r 111206_0150_BD0BP3ACXX_1.locus.pileup.rpt -p 111206_0150_BD0BP3ACXX_1.locus.pileup.prop
```

Calculate RAD locus summary statistics:
```
summarise_rad_loci.pl -t 111206_0150_BD0BP3ACXX_1.locus.pileup.csv -r 111206_0150_BD0BP3ACXX_1.locus.pileup.rpt -p 111206_0150_BD0BP3ACXX_1.locus.pileup.prop -f c_elegans.WS229.genomic.PstI.fragments.txt -l 96 > 111206_0150_BD0BP3ACXX_1.locus_summary.tsv
```

These output files could be used for many other analyses, and many statistics are unused in the paper. The `summarise_rad_loci.pl` script will only work correctly with output from paired end data but `rad_locus_perfect_matches.pl` should work with single or paired end data.


The plot is generated as follows:

```r
> library(ggplot2)
> library(scales)
> N2.uniqtags<-read.delim("111206_0150_BD0BP3ACXX_1.locus_summary.tsv")
> N2.uniqtags.14.1<-N2.uniqtags[N2.uniqtags$PCR==14 & N2.uniqtags$Rep==1,]
> pdf("Figure3.pdf")
> ggplot(N2.uniqtags.14.1[N2.uniqtags.14.1$SampleReads>0,], aes(GenomeFragLen,SampleReads))+geom_point(size=0.5)+theme_bw()+scale_x_log10(name="Restriction Fragment Length",breaks=c(1,10,100,1000,10000,100000),labels=comma,limits=c(100,100000))+scale_y_continuous(name="Read Depth",limits=c(0,700))+theme_bw(base_family="Palatino")
> dev.off()
```

### Correlations

The correlations reported in the Results section "Restriction fragment length biases read depth at RAD loci" were calculated using the TSV file generated above as follows:

````r
> N2.uniqtags<-read.delim("111206_0150_BD0BP3ACXX_1.locus_summary.tsv")
> N2.uniqtags.14.1<-N2.uniqtags[N2.uniqtags$PCR==14 & N2.uniqtags$Rep==1,]
> cor.test(log(N2.uniqtags.14.1[N2.uniqtags.14.1$SampleReads>0,]$GenomeFragLen),N2.uniqtags.14.1[N2.uniqtags.14.1$SampleReads>0,]$SampleReads)

	Pearson's product-moment correlation

data:  log(N2.uniqtags.14.1[N2.uniqtags.14.1$SampleReads > 0, ]$GenomeFragLen) and N2.uniqtags.14.1[N2.uniqtags.14.1$SampleReads > 0, ]$SampleReads 
t = 161.7588, df = 24826, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0 
95 percent confidence interval:
 0.7102265 0.7223394 
sample estimates:
      cor 
0.7163369 

> cor.test(log(N2.uniqtags.14.1[N2.uniqtags.14.1$SampleReads>0 & N2.uniqtags.14.1$GenomeFragLen>10000,]$GenomeFragLen),N2.uniqtags.14.1[N2.uniqtags.14.1$SampleReads>0 & N2.uniqtags.14.1$GenomeFragLen>10000,]$SampleReads)

	Pearson's product-moment correlation

data:  log(N2.uniqtags.14.1[N2.uniqtags.14.1$SampleReads > 0 & N2.uniqtags.14.1$GenomeFragLen >  and N2.uniqtags.14.1[N2.uniqtags.14.1$SampleReads > 0 & N2.uniqtags.14.1$GenomeFragLen >      10000, ]$GenomeFragLen) and     10000, ]$SampleReads 
t = 3.0396, df = 6774, p-value = 0.002378
alternative hypothesis: true correlation is not equal to 0 
95 percent confidence interval:
 0.01310666 0.06066359 
sample estimates:
       cor 
0.03690602 

> cor.test(log(N2.uniqtags.14.1[N2.uniqtags.14.1$SampleReads>0 & N2.uniqtags.14.1$GenomeFragLen<=10000,]$GenomeFragLen),N2.uniqtags.14.1[N2.uniqtags.14.1$SampleReads>0 & N2.uniqtags.14.1$GenomeFragLen<=10000,]$SampleReads)

	Pearson's product-moment correlation

data:  log(N2.uniqtags.14.1[N2.uniqtags.14.1$SampleReads > 0 & N2.uniqtags.14.1$GenomeFragLen <=  and N2.uniqtags.14.1[N2.uniqtags.14.1$SampleReads > 0 & N2.uniqtags.14.1$GenomeFragLen <=      10000, ]$GenomeFragLen) and     10000, ]$SampleReads 
t = 122.4107, df = 18050, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0 
95 percent confidence interval:
 0.6654481 0.6813912 
sample estimates:
     cor 
0.673498
````


### Figure 4

Using RAD locus summary file output by `summarise_rad_loci.pl` above:

```r
> N2.uniqtags<-read.delim("111206_0150_BD0BP3ACXX_1.locus_summary.tsv")
> N2.uniqtags.gt10K<-N2.uniqtags[N2.uniqtags$GenomeFragLen>=10000,]
> pdf("Figure4A.pdf",10,2)
> ggplot(N2.uniqtags.gt10K[N2.uniqtags.gt10K$Rep==2 & N2.uniqtags.gt10K$SampleReadsNorm>0,],aes(MeanTagGC,SampleReadsNorm))+facet_grid(.~PCR)+geom_point(size=0.5)+ylim(c(0,100))+geom_smooth(method="loess")+theme_bw(base_family="Palatino",base_size=10)+ylab("Normalized Read Depth")+xlab("Mean Sheared Fragment GC Content (%)")
> dev.off()
> pdf("Figure4B.pdf",10,2)
> ggplot(N2.uniqtags.gt10K[N2.uniqtags.gt10K$Rep==2 & N2.uniqtags.gt10K$SampleFragmentsNorm>0,],aes(MeanTagGC,SampleFragmentsNorm))+facet_grid(.~PCR)+geom_point(size=0.5)+ylim(c(0,100))+geom_smooth(method="loess")+theme_bw(base_family="Palatino",base_size=10)+ylab("Normalized Fragment Depth")+xlab("Mean Sheared Fragment GC Content (%)")
> dev.off()
```



### Figure S1

This figure uses public stickleback SbfI RAD run [SRR034312](http://www.ncbi.nlm.nih.gov/sra?term=SRR034312) published in Hohenlohe et al., PLoS Genetics, 2010, and C. elegans EcoRI RAD library ECA17, Sequence Read Archive accession 047839 (submitted but not public at time of writing, kindly supplied by Erik Andersen and Josh Shapiro), published in Andersen et al., Genetics, 2012.

To get stickleback FASTQ data, use the `fastq-dump` tool from the [SRA Toolkit](http://www.ncbi.nlm.nih.gov/books/NBK47540/):

````
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR034/SRR034312/SRR034312.sra
fastq-dump SRR034312.sra
````

Simulate RAD fragments for the stickleback genome ([assembly BROAD S1 in Ensembl release 66](ftp://ftp.ensembl.org/pub/release-66/fasta/gasterosteus_aculeatus/dna/)):

````
simulate_rad_fragments.pl -i Gasterosteus_aculeatus.BROADS1.66.dna.chromosomes.fa -e CCTGCAGG > Gasterosteus_aculeatus.BROADS1.66.dna.chromosomes.SbfI.fragments.txt &
````

Simulate EcoRI fragments for the C. elegans genome:

````
simulate_rad_fragments.pl -i c_elegans.WS229.genomic.fa -e GAATTC > c_elegans.WS229.EcoRI.fragments.txt
````


Generate RAD locus CSV files for each library, truncating the two C. elegans libraries to match the number of passed reads from the stickleback (the smallest library):

````
rad_locus_perfect_matches.pl -f c_elegans.WS229.genomic.PstI.fragments.txt -l 96 -e C*TGCA*G -1 <(gunzip -c 111206_0150_BD0BP3ACXX_1_1.sanfastq.gz) -b celegans_rad_samples.pools -m 6516557 -c 111206_0150_BD0BP3ACXX_1.trunc.tag.pileup.uniqlocs.all.csv -r 111206_0150_BD0BP3ACXX_1.trunc.tag.pileup.uniqlocs.all.rpt -p 111206_0150_BD0BP3ACXX_1.trunc.tag.pileup.uniqlocs.all.prop

rad_locus_perfect_matches.pl -f c_elegans.WS229.EcoRI.fragments.txt -l 96 -e G*AATT*C -1 <(cat s_2_sequence.txt s_3_sequence.txt) -b N2-like.pools -m 6516557 -c N2-like.trunc.tag.pileup.uniqlocs.all.csv -r N2-like.trunc.tag.pileup.uniqlocs.all.rpt -p N2-like.trunc.tag.pileup.uniqlocs.all.prop 

rad_locus_perfect_matches.pl -f Gasterosteus_aculeatus.BROADS1.66.dna.chromosomes.SbfI.fragments.txt -l 49 -e CC*TGCA*GG -1 SRR034312.fastq -b SRR034312.pools -c SRR034312.tag.pileup.uniqlocs.all.csv -r SRR034312.tag.pileup.uniqlocs.all.rpt -p SRR034312.tag.pileup.uniqlocs.all.prop
````


Make plots for each CSV file, for example, for stickleback:

````r
> read.csv("SRR034312.tag.pileup.uniqlocs.all.csv")->stickleback.fragstats
> pdf("Supp1_stickleback.FragLen.TotalDepth.pdf",7,5)
> ggplot(stickleback.fragstats[stickleback.fragstats$TotalDepth>0,], aes(ResFragLen,TotalDepth))+geom_point(size=0.5)+theme_bw()+scale_x_log10(name="Restriction Fragment Length",breaks=c(300,500,700,1000,10000,100000,200000),labels=comma,limits=c(300,200000))+scale_y_continuous(name="Read Depth",limits=c(0,600))+theme_bw(base_family="Palatino")
> dev.off()

> pdf("Supp1_stickleback.SbfI.FragLen.histogram.pdf",7,5)
> ggplot(stickleback.fragstats,aes(ResFragLen))+geom_histogram(binwidth=0.05)+scale_x_log10(name="Restriction Fragment Length",breaks=c(300,500,700,1000,10000,100000,200000),labels=comma,limits=c(300,200000))+scale_y_continuous(name="Frequency",limits=c(0,3600))+theme_bw(base_family="Palatino")
> dev.off()
````


### Figure S3

Generate RAD locus CSV files for the two libraries N2rad14 ("Single") and N2rad15 ("Mixed"), with N2rad15 truncated to match the passed read numbers from N2rad14:

````
rad_locus_perfect_matches.pl -f /ifs/blaxterprojects/2011151_Davey_John/c_elegans.WS229.genomic.PstI.fragments.txt -l 143 -e C*TGCA*G -1 <(cat N2rad14*_1.fastq) -b N2rad14.pools -c N2rad14.tag.pileup.uniqlocs.all.csv -r N2rad14.tag.pileup.uniqlocs.all.rpt -p N2rad14.tag.pileup.uniqlocs.all.prop
rad_locus_perfect_matches.pl -f /ifs/blaxterprojects/2011151_Davey_John/c_elegans.WS229.genomic.PstI.fragments.txt -l 143 -e C*TGCA*G -1 <(cat N2rad15*_1.fastq) -b N2rad15.pools -c N2rad15.tag.pileup.uniqlocs.all.csv -r N2rad15.tag.pileup.uniqlocs.all.rpt -p N2rad15.tag.pileup.uniqlocs.all.prop -m 1396418
````

Plot figure using the two CSV files:

````r
> library(ggplot2)
> library(scales)
> read.csv("N2rad14.tag.pileup.uniqlocs.all.csv")->N2rad14.fragstats
> read.csv("N2rad15.tag.pileup.uniqlocs.all.csv")->N2rad15.fragstats
> N2rad14.join<-data.frame(ResFragLen=N2rad14.fragstats$ResFragLen,TotalDepth=N2rad14.fragstats$TotalDepth,Library="Mixed")
> N2rad15.join<-data.frame(ResFragLen=N2rad15.fragstats$ResFragLen,TotalDepth=N2rad15.fragstats$TotalDepth,Library="Single")
> N2rad.join<-rbind(N2rad14.join,N2rad15.join)
> pdf("N2rad14+15.fragstats.pdf",10,7)
> ggplot(N2rad.join[N2rad.join$TotalDepth>0,], aes(ResFragLen,TotalDepth,colour=factor(Library)))+geom_point(size=0.25)+theme_bw()+scale_x_log10(name="Restriction Fragment Length",breaks=c(100,1000,10000,50000),labels=comma,limits=c(100,50000))+scale_y_continuous(name="Read Depth",limits=c(0,150))+theme_bw(base_family="Palatino")+geom_smooth(method="loess")+scale_colour_hue(name="Shearing")
> dev.off()
````



## HELICONIUS DATA SET ##

For the Heliconius data, a pipeline was developed to produce a similar summary TSV file as that used for Figure 3. However, this file was produced by inferring RAD loci from an alignment, rather than identifying them in the genome first. Also, haplotypes for the aligned RAD cross were inferred to assess allele copy number. This pipeline has not been tested on any other data and will only work with paired end RAD data.

FASTQs from Heliconius melpomene melpomene x Heliconius melpomene rosina individuals from ENA project [ERP000993](http://www.ebi.ac.uk/ena/data/view/ERP000993) and [Heliconius melpomene reference genome v1.1](http://metazoa.ensembl.org/Heliconius_melpomene) were used for this analysis.

### Alignment

Prepare genome for use with Stampy: 

```
stampy.py -G Hmel1-1_primaryScaffolds Hmel1-1_primaryScaffolds.fa
stampy.py -H Hmel1-1_primaryScaffolds -g Hmel1-1_primaryScaffolds
```

Align each offspring separately using Stampy, parallelised using an SGE pipeline on the Edinburgh compute cluster [Eddie](http://www.ecdf.ed.ac.uk/). This pipeline is highly specific to the Edinburgh cluster and so is not publicly available. Please contact us for more details if interested. For the purposes of the RAD analysis, we believe any BAM file from any aligner should be fine, provided it conforms to the SAM specification. But as the downstream scripts process CIGAR fields among other things, it would not be surprising to find other aligners output slightly different formats that confuse the scripts.

The command below shows typical Stampy options for one offspring, in this case number 1. The total number of parts changes depending on the number of reads sequenced for each individual.

```
stampy.py --gatkcigarworkaround --baq  --alignquals --insertsize=500 --insertsd=100 --processpart $SGE_TASK_ID/24 --readgroup=ID:PstI.1.110802,SM:PstI.1,PL:Illumina -g ./Hmel1-1_primaryScaffolds -h ./Hmel1-1_primaryScaffolds -o PstI.1.110802.Hmel1-1_primaryScaffolds.$SGE_TASK_ID.sam -M 1_1.fastq 1_2.fastq
java -Xmx1500m -jar $PICARDPATH/SortSam.jar TMP_DIR=$TMPDIR MAX_RECORDS_IN_RAM=50000 SORT_ORDER=coordinate INPUT=PstI.1.110802.Hmel1-1_primaryScaffolds.$SGE_TASK_ID.sam OUTPUT=PstI.1.110802.Hmel1-1_primaryScaffolds.$SGE_TASK_ID.bam
```

Merge all part BAM files for a particular individual, adding an INPUT argument for each BAM file generated by parts:

```
java -Xmx7g -jar $PICARDPATH/MergeSamFiles.jar TMP_DIR=$TMPDIR CREATE_INDEX=true MAX_RECORDS_IN_RAM=50000 OUTPUT=PstI.1.110802.Hmel1-1_primaryScaffolds.bam INPUT=PstI.1.110802.Hmel1-1_primaryScaffolds.1.bam INPUT=PstI.1.110802.Hmel1-1_primaryScaffolds.2.bam INPUT=....
```


Then merge all BAM files for all individuals into one master BAM file:

```
java -Xmx7g -jar $PICARDPATH/MergeSamFiles.jar CREATE_INDEX=true OUTPUT=PstI.all.110802.Hmel1-1_primaryScaffolds.bam INPUT=PstI.1.110802.Hmel1-1_primaryScaffolds.bam INPUT=PstI.2.110802.Hmel1-1_primaryScaffolds.bam INPUT=PstI.3.110802.Hmel1-1_primaryScaffolds.bam INPUT=....
```

The BAM file for just the used chromosome 18 scaffolds is provided as `PstI.all.110802.Hmel1-1_primaryScaffolds.chr18.bam` in European Nucleotide Archive project ERP001757, experiment ERX144563.


### Infer RAD loci

Infer RAD loci reads per individual from the BAM file. This script outputs a different CSV file for each scaffold in the genome to a folder named as the output parameter. It is very inefficient and takes ~1 day to run on the Heliconius BAM. The CSV output files are analogous to those produced by `rad_locus_perfect_matches.pl` described above.

````
infer_rad_loci_per_scaffold_from_bam.pl -b PstI.all.110802.Hmel1-1_primaryScaffolds.chr18.bam -o stampy.csv
````

Then generate TSV summary files from the CSV files, passing tag length (ie read length after barcode is removed) as an argument. For the Heliconius data this is 95.

````
for i in stampy.csv/*csv; do (infer_rad_loci_statistics.pl -l $i -r Hmel1-1_primaryScaffolds.fa -t 95 -o ${i%csv}tsv); done
````

This produces TSV files similar to those produced for C. elegans.


### Infer haplotypes

The previous scripts are intended to be general enough to work on other data sets (which is not to say that they will work). These haplotyping scripts, however, are tightly integrated with the Heliconius cross and are not intended for use on other data sets. They depend on a set of genotypes and other cleaning information that was manually generated for this cross; they also have hardcoded information about the parents and the cross, which are almost certainly going to be different for another cross. This script is provided for transparency and direct replication only.

For each scaffold of interest:

````
infer_rad_loci_haplotypes.pl -s scaffold_name
````

The `infer_rad_loci_haplotypes.pl` script will expect the following files to be available in the working directory for scaffold `scaffold_name`:

PstI.all.110802.Hmel1-1_primaryScaffolds.chr18.bam.scaffold_name.tsv (produced by infer_rad_loci_statistics.pl)
genotypes.scaffold_name.tsv
tippex.scaffold_name.tsv

More details on the `genotypes` and `tippex` formats can be found below. Files for the Heliconius cross and loci used are included in this distribution. 


The `plot_restriction_sites_fraglen.R` script can be used to visualize the RAD loci aligned to the scaffold for all individuals. This visualization was used to identify valid RAD loci and ignore restriction fragments from loci that could not be inferred from the alignment.

````
plot_restriction_sites_fraglen.R scaffold_name
````

To run for usable chromosome 18 Heliconius scaffolds:

````
perl -ne 'chomp;my ($scf,$length)=split /\t/; system("infer_haplotypes_from_inferred_rad_loci_recombinants.pl -s $scf");system("../plot_restriction_sites_fraglen.R $scf")' chr18.scaffolds.clean.names.tsv
````

The plots created by this command are in the heliconius_scaffold_files folder (filenames scaffold_name.ressites.fragments.pdf).

Merge output for these scaffolds:

````
perl -ne 'BEGIN{$first = 1}{chomp; open $haps, "<", "PstI.all.110802.Hmel1-1_primaryScaffolds.chr18.bam.$_\.haplotypes.tsv"; $header=<$haps>; if ($first) {print $header; $first=0}; print <$haps>; close $haps;}' chr18.scaffolds.clean.names.tsv > chr18.scaffolds.clean.loci.tsv
````

### Figure 6

Figure 6 is generated using the `chr18.scaffolds.clean.loci.tsv` file generated above and the following commands:

````
figure_6_data.pl -t chr18.scaffolds.clean.loci.tsv > chr18.scaffolds.clean.missing.data.tsv
figure_6_plot.R chr18.scaffolds.clean.missing.data.tsv psti.rad.read.counts.tsv chr18.scaffolds.clean.loci.tsv
````

The Perl script generates summary statistics about loci with missing genotypes. The R script plots all three parts of the figure on one plot, output to a PDF file. The `psti.rad.read.counts.tsv` file is a simple list of read counts per individual in the Heliconius PstI cross (included in this distribution).

### Figure S4 

This figure is a modification of the top plot of Figure 6 and requires the `chr18.scaffolds.clean.loci.tsv` file generated for that figure.

````r
> library(ggplot2)
> library(scales)
> chr18.clean<-read.delim("chr18.scaffolds.clean.loci.tsv")
> samples<-c("F1Mo","F1Fa","28","24","73","75")
> chr18.clean.select<-chr18.clean[chr18.clean$Sample %in% samples & chr18.clean$AlleleCount>0 & (!is.na(chr18.clean$MoFragLen) | !is.na(chr18.clean$FaFragLen) | !is.na(chr18.clean$SharedFragLen)),]
> chr18.clean.select<-transform(chr18.clean.select,Sample=factor(chr18.clean.select$Sample,levels=c("F1Mo","F1Fa","28","24","73","75"), labels=c("Mother (AC)", "Father (BC)", "AB","AC","CB","CC")))
> pdf("figure_S4_heliconius_fragment_depth.pdf",13.30,3.6)
> ggplot(chr18.clean.select,
        aes(SharedFragLen,SampleAlleleFragments,colour=factor(AlleleCount)))+
        facet_grid(.~Sample)+
        geom_point(size=1)+
        ylim(c(0,100))+
        xlim(c(1,1200000))+
        scale_x_continuous(trans=log10_trans(),breaks=c(400,1000,4000,10000,40000),labels=c("0.4","1","4","10","40"))+
        theme_bw(base_family="Palatino")+
        scale_colour_manual(name="Allele\nCopies",values=c("mediumorchid1","forestgreen"))+
        xlab("Restriction Fragment Length (kbp)")+
        ylab("Sheared Fragment Depth")
> dev.off()
````


### Figure 7

Tool comparison requires selection of reads that contribute to curated haplotypes, listed in `haplotypes.tsv` files generated above from the master BAM file `PstI.all.110802.Hmel1-1_primaryScaffolds.chr18.bam`.

Extract reads from cleaned scaffolds into BAM files:

````
mkdir toolcomp
perl -ne 'chomp; system("samtools view -bo toolcomp/PstI.all.110802.Hmel1-1_primaryScaffolds.$_\.bam PstI.all.110802.Hmel1-1_primaryScaffolds.chr18.bam $_")' chr18.scaffolds.clean.names.tsv
````

Generate SAM files for each cleaned scaffold and FASTQ files for each individual, for processing by Stacks, RADtools and GATK. SAM files will be written to `toolcomp/toolcomp.samfiles` and FASTQ files to `toolcomp/toolcomp.fastqfiles`.

````
perl -ne 'chomp; system("prepare_rad_tool_input_files.pl -b toolcomp/PstI.all.110802.Hmel1-1_primaryScaffolds.$_\.bam -t PstI.all.110802.Hmel1-1_primaryScaffolds.chr18.bam.$_\.haplotypes.tsv -s $_ -o toolcomp");' chr18.scaffolds.clean.names.tsv 2> chr18.scaffolds.clean.extract.reads.log
````

Merge SAM files, adding header from master BAM file, sort, then select read 1 only so GATK results can be compared to RADtools and Stacks (as the de novo tools only cluster on read 1 and do not use read 2 information):

````
cat <(samtools view -H PstI.all.110802.Hmel1-1_primaryScaffolds.chr18.bam) toolcomp/toolcomp.samfiles/*sam > toolcomp/PstI.all.Hmel1-1_primaryScaffolds.chr18.scaffolds.clean.unsorted.sam
java -jar $PICARDPATH/SortSam.jar INPUT=toolcomp/PstI.all.Hmel1-1_primaryScaffolds.chr18.scaffolds.clean.unsorted.sam OUTPUT=toolcomp/PstI.all.Hmel1-1_primaryScaffolds.chr18.scaffolds.clean.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true
samtools view -b -f 64 -o toolcomp/PstI.all.Hmel1-1_primaryScaffolds.chr18.scaffolds.clean.sorted.read1.bam toolcomp/PstI.all.Hmel1-1_primaryScaffolds.chr18.scaffolds.clean.sorted.bam
samtools index toolcomp/PstI.all.Hmel1-1_primaryScaffolds.chr18.scaffolds.clean.sorted.read1.bam
````

Running tools:

#### GATK

````
java -jar $GATK_PATH/GenomeAnalysisTK.jar -T UnifiedGenotyper -I toolcomp/PstI.all.Hmel1-1_primaryScaffolds.chr18.scaffolds.clean.sorted.read1.bam -R Hmel1-1_primaryScaffolds.fa -o toolcomp/PstI.all.Hmel1-1_primaryScaffolds.chr18.scaffolds.clean.sorted.read1.GATK.vcf -out_mode EMIT_ALL_CONFIDENT_SITES -glm BOTH
````

#### RADtools

(Select an appropriate `max_processes` value depending on the number of cores on your machine.)

````
cp toolcomp.fastqfiles.pools toolcomp/toolcomp.fastqfiles
cd toolcomp/toolcomp.fastqfiles
[toolcomp/toolcomp.fastqfiles] $ RADtags -z -s TGCAG --max_processes 47 -v 2>&1 >> RADtags.log &
[toolcomp/toolcomp.fastqfiles] $ for i in toolcomp*; do mv $i ${i#toolcomp.}; done
[toolcomp/toolcomp.fastqfiles] $ RADmarkers -v > RADmarkers.output 2> RADmarkers_revision.log &
````


#### Stacks

Stacks could be run with the de novo mapping pipeline or a pipeline of your choice, but the component command lines are given here to show the exact options used.


Set up directory structure for Stacks and create links to input FASTQ files:
````
[toolcomp] mkdir stacks
[toolcomp] mkdir stacks/samples
[toolcomp] mkdir stacks/stacks
[toolcomp/stacks/samples] for i in ../../toolcomp.revision.fastqfiles/*_1.fastq; do j=`echo $i | awk 'BEGIN{FS="/"}{print $NF}'`; ln -s $i $j ; done
````

For each individual, run `ustacks` (eg for the F1 father F1Fa) with 5 mismatches (to match RADtools defaults):

````
[toolcomp/stacks] ustacks -t fastq -f F1Fa_1.fastq -o stacks -i 1 -d -r -M 5 -m 2
````

Build a catalog:

````
[toolcomp] cstacks -b 1 -o stacks -s stacks/F1Fa.tags.tsv -S 1 -s stacks/F1Mo.tags.tsv -S 2 -n 2
````

Run sstacks against the catalog for each individual, eg for F1Fa:

````
[toolcomp] sstacks -b 1 -c stacks/batch_1 -s F1Fa -S 1 -o stacks
````

Generate genotypes with corrections:

````
[toolcomp] genotypes -P stacks -b 1 -t gen -c -r 1
````

#### Plot

The `figure_7_data.pl` script compiles and validates the output of GATK, RADtools and Stacks, and the `figure_7_plot.R` script generates a plot based on this data:

````
figure_7_data.pl -t chr18.scaffolds.clean.loci.tsv -r toolcomp/toolcomp.fastqfiles -v toolcomp/PstI.all.Hmel1-1_primaryScaffolds.chr18.scaffolds.clean.sorted.read1.GATK.vcf -m toolcomp/toolcomp.fastqfiles/RADmarkers.output -s toolcomp/stacks/stacks -c toolcomp/stacks/stacks/batch_1.catalog.tags.tsv -h toolcomp/stacks/stacks/batch_1.haplotypes_1.tsv  > figure_7.tools.txt 2> figure_7.loci.txt
figure_7_plot.R figure_7.tools.txt figure_7_readprops.pdf figure_7_fragment_skew.pdf
````


## FILE FORMATS


### pools

These files have two columns, separated by any whitespace, with sample name in the first column and barcode in the second. Multiple barcodes can be given for each sample. See RADtools documentation for more information.

### psti.rad.read.counts.tsv

This file lists Sample names and Reads aligned to the whole genome for the Heliconius PstI cross in two columns.

### chr18.scaffolds.clean.names.tsv

This file is a simple list of the chromosome 18 Heliconius scaffolds containing usable RAD loci.

### fragments

The `c_elegans.WS229.genomic.PstI.fragments.txt` file is included for reference with regard to the `c_elegans.loci.tsv` file, as the in silico RAD loci are derived from the in silico RAD fragments. The file contains one line for each restriction fragment, described by six columns:

````
Chromosome - name of chromosome featuring restriction fragment
Start      - genomic start position of restriction fragment
End        - genomic end position of restriction fragment
Length     - length of restriction fragment
GC         - GC content of restriction fragment
Sequence   - full sequence of restriction fragment
````

All fragments retain full restriction sites at the start and end, except for those at the start and end of each chromosome.

### loci

The two loci files `c_elegans.loci.tsv` and `chr18.scaffolds.clean.loci.tsv` are related but distinct formats for in silico and inferred loci. Each line reports summary statistics for one sample aligned to one RAD locus. Some fields report general statistics for the RAD locus over all samples.

#### c_elegans.loci.tsv ####

````
Chr                 - chromosome of RAD locus
Loc                 - position of RAD locus on chromosome
GenomeFragLen       - length of restriction fragment containing RAD locus
MinTagFragLen       - shortest sheared fragment length aligned to RAD locus for all samples
MaxTagFragLen       - longest sheared fragment length aligned to RAD locus for all samples
MinTagLenGC         - GC content of shortest sheared fragment aligned to RAD locus for all samples
MeanTagGC           - mean GC content of all sheared fragments aligned to RAD locus for all samples
MaxTagLenGC         - GC content of longest sheared fragment aligned to RAD locus for all samples
TotalReads          - all reads aligned to this RAD locus from all samples
Sample              - sample name
PCR                 - PCR cycles used for sample (for C. elegans library, otherwise NA)
Rep                 - Replicate used for sample (for C. elegans library, otherwise NA)
SampleFragments     - unique fragments (after PCR duplicates removed) aligned to RAD locus for this sample
SampleFragmentsNorm - normalised fragment count (see Methods for calculation)
SampleReads         - reads aligned to RAD locus for this sample
SampleReadsNorm     - normalised read count (see Methods for calculation)
MinSampleFragLen    - shortest sheared fragment length aligned to RAD locus for this sample
MeanSampleFragLen   - mean sheared fragment length aligned to RAD locus for this sample
MaxSampleFragLen    - longest sheared fragment length aligned to RAD locus for this sample
MinSampleLenGC      - GC content of shortest sheared fragment aligned to RAD locus for this sample
MeanSampleLenGC     - mean GC content of all sheared fragments aligned to RAD locus for this sample
MaxSampleLenGC      - GC content of longest sheared fragment aligned to RAD locus for this sample
````

#### chr18.scaffolds.clean.loci.tsv ####

The Heliconius loci have many additional fields, due to two major features - the loci are heterozygous, and therefore multiple alleles are reported for each locus; and haplotypes have been called for each locus based on the Heliconius mapping cross. These fields have the same meaning as above:

````
Chr
Loc
MinTagFragLen
MaxTagFragLen
MinTagLenGC
MaxTagLenGC
Sample
MinSampleFragLen
MeanSampleFragLen
MaxSampleFragLen
MinSampleLenGC
MeanSampleLenGC
MaxSampleLenGC
````

These fields are new:

````
Dir                   - Direction of the RAD locus (+/-)
Allele                - number of alleles (1,2,3 can be genuine; 4 and above are usually repeats or sequencing errors)
Seq                   - sequence of allele
TotalLocusReads       - total reads for the RAD locus over all samples
TotalAlleleReads      - total reads for the allele over all samples
SampleLocusReads      - reads for the RAD locus for this sample
SampleLocusFragments  - fragments (PCR duplicates removed) for the RAD locus for this sample
SampleAlleleReads     - reads for the allele for this sample
SampleAlleleFragments - fragments for the allele for this sample
SampleAlleleReadsProp - proportion of locus reads contributing to this allele (towards 1 usually indicates a homozygous allele, towards 0.5 heterozygous)
SampleAlleleFragsProp - proportion of locus fragments contributing to this allele
SampleResFragLen      - DO NOT USE - crude estimate of fragment length, superseded by fields below
SampleMatchingTag     - DO NOT USE - crude estimate of neighbouring RAD locus, superseded by fields below
AlleleCount           - copies of this allele present in sample based on haplotyping
MoHap                 - maternal haplotype (F1MoA or F1MoC for offspring)
MoFragLen             - fragment length of maternal haplotype, given neighbouring RAD locus
MoSite                - position of maternal RAD locus (same as Loc)
MoNeighbour           - neighbouring maternal RAD locus (next locus in direction Dir, ie at other end of restriction fragment)
FaHap                 - paternal haplotype (F1FaB or F1FaC for offspring)
FaFragLen             - fragment length of paternal haplotype, given neighbouring RAD locus
FaSite                - position of paternal RAD locus (same as Loc)
FaNeighbour           - neighbouring paternal RAD locus (next locus in direction Dir, ie at other end of restriction fragment)
SharedFragLen         - fragment length when MoFragLen and FaFragLen are equal, otherwise NA
````

### genotypes

The `genotypes` files, one for each scaffold examined, define the master genotypes along each scaffold for each offspring in the Heliconius RAD cross. Parental genotypes are defined in terms of the three primary haplotypes A, B and C:

````
Sample  Type    Start   End     Mother  Father
F1Mo    P       1       1325285 A       C
F1Fa    P       1       1325285 B       C
````

Offspring genotypes are then defined in relation to the parental genotypes, where 1 is the Mother haplotype (A or B in this case) and 2 is the Father haplotype (both C here).

````
F0GM    GP      1       1325285 1       1
1       O       1       1325285 1       1
6       O       1       1325285 1       2
24      O       1       1325285 1       2
19      O       1       1325285 2       1
33      O       1       1325285 2       1
109     O       1       1325285 2       2
116     O       1       1325285 2       2
````

Here, for example, offspring 24 has genotype AC and offspring 33 has genotype CB.

The type field is either P (parent), GP (grandparent) or O (offspring).

If an individual has a recombination along the scaffold, this is specified using two lines and defining appropriate start and end points:

````
26      O       1       952000  2       1
26      O       952001  1325285 2       2
58      O       1       325000  2       1
58      O       325001  1325285 2       2
````

These points were defined manually based on the plots generated by the `plot_restriction_sites_fraglen.R` script.


### tippex

The `tippex` files define the sites to be ignored or corrected because of low coverage or repeat content. Four field types are defined

````
min	1
max	100000
````

These values define the region of a scaffold with acceptable restriction fragments. Often the ends of scaffolds are repeat content (which is why the scaffold can not be continued), so many sites at the ends of scaffolds need to be ignored.

````
ignore	50000
````

This causes downstream scripts to ignore the RAD locus at 50000 bp.

````
hap     87051   +       F1FaB   3
````

This specifies a missing haplotype, which has not been inferred automatically because at least one individual has no coverage for this allele. The second and third fields define the position and direction of the RAD locus on the scaffold; the fourth field is the haplotype that needs correcting, and the fifth field defines the allele at the locus that should be assigned to the haplotype.

Fields can be commented out with a hash symbol.


### ressites.fragments.pdf ###

A plot is provided for each used chromosome 18 Heliconius scaffold. These plots were used as a visual aid to curate the usable RAD loci and identify recombinations, although only the final used RAD loci are present on the current PDFs. Re-run the commands without supplying tippex files to see the whole scaffold.

Scaffold position is shown from left to right. Each row shows coverage for a particular individual, parents at the top, offspring below. For each row, the band at the base names the individual and genotype, with genotypes coloured as follows:

````
AC  Red   F1 Mother
BC  Blue  F1 Father
AB  Green F0 Grandmother

````
