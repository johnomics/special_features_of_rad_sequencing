#!/usr/bin/env Rscript

# plot_restriction_sites.R
# John W Davey john.davey@ed.ac.uk
# Load Heliconius colour pattern restriction site data from TSV file
# and output positions and depths of sites in all individuals

plot.sample<-function(sample,sample.sites.all,samples,genotypes,scaffold.length, scaffold.agp) {
    
    num.samples<- length(samples)
    sample.frac <- 1/num.samples
    sample.ind <- which(samples==sample)

    sample.sites<-sample.sites.all[sample.sites.all$Allele<=4,]

    # Get Brewer colours for each genotype
    genotypes$Genotype<-paste(genotypes$Mother,genotypes$Father,sep="")
    
    gt.col<-data.frame(Genotype=unique(genotypes$Genotype),Colour=brewer.pal(length(unique(genotypes$Genotype)),"Set1"))
    
    # Sample viewport
    sample.vp<-viewport(0,sample.ind/num.samples,just=c("left","bottom"),width=1,height=1/num.samples)
    pushViewport(sample.vp)
    
    # Label viewport
    
    label.vp<-viewport(0,0,just=c("left","bottom"),width=0.1,height=1)
    pushViewport(label.vp)
    grid.text(sample,just="left")
    popViewport() #label.vp
    
    # Scaffold viewport
    scaffold.vp<-viewport(0.1,0,just=c("left","bottom"),width=0.9,height=1,xscale=c(0,scaffold.length),yscale=c(0,8))
    pushViewport(scaffold.vp)

    # Scaffold

    # Get current genotype and plot stripe labelling genotype for individual by name and by colour
    cur.gt<-genotypes[genotypes$Sample==sample,]
    grid.polyline(
        unit(c(cur.gt$Start,cur.gt$End),"native"),
        unit(rep(c(0,0),length(cur.gt$Start)),"native"),
        id=rep(1:length(cur.gt$Start),2),
        gp=gpar(lwd=1,lineend="round",col=rep(sapply(cur.gt$Genotype,function(x){as.character(gt.col[gt.col$Genotype==x,]$Colour)}),2))
    )
    # Plot individual name along stripe
    label.pos<-seq(1000,scaffold.length,by=20000)
    grid.text(sample,unit(label.pos+10000,"native"),unit(0,"native"),gp=gpar(fontsize=1,col="black"))

    # Plot genotypes (including recombinants)
    sapply(
        cur.gt$Genotype,
        function(x) {
            grid.text(
                x,
                unit(label.pos[label.pos > cur.gt[cur.gt$Genotype==x,]$Start & label.pos < cur.gt[cur.gt$Genotype==x,]$End],"native"),
                unit(0,"native"),
                gp=gpar(fontsize=1,col="black")
            )
        }
    )
    
    num.alleles<-length(sample.sites$Loc)

    loci.dir<-sapply(sample.sites$Dir,function(x) {if (x=="+") 1 else -1})

    max.depth<-max(sample.sites$SampleAlleleFragments)

    # Loci with arrows
    grid.polyline(
        unit(c(sample.sites$Loc,sample.sites$Loc+(1000*loci.dir)),"native"),
        rep(unit(3.5+sample.sites$Allele,"native"),2),
        arrow=arrow(length=unit(0.02,"inches"),type="closed"),
#        arrow=arrow(length=unit(500,"native"),type="closed"),
        id=rep(1:num.alleles,2),
        gp=gpar(lwd=0.75, alpha=rep(sample.sites$SampleAlleleFragments/max.depth,2),col=rep(sample.sites$Allele,2))
    )

    # Read 1
    grid.polyline(
        unit(c(sample.sites$Loc,sample.sites$Loc+(96*loci.dir)),"native"),
        rep(unit(3.45+sample.sites$Allele,"native"),2),
        id=rep(1:num.alleles,2),
        gp=gpar(lwd=0.25,col=rep(sample.sites$Allele,2))
    )

    # Read 2 contig
    grid.polyline(
        unit(c(sample.sites$Loc+(sample.sites$MinSampleFragLen*loci.dir),sample.sites$Loc+(sample.sites$MaxSampleFragLen*loci.dir)),"native"),
        rep(unit(3.55+sample.sites$Allele,"native"),2),
        id=rep(1:num.alleles,2),
        gp=gpar(lwd=0.25,col=rep(sample.sites$Allele,2))
    )

    # Loci positions in text
    loci.text.just<-sapply(sample.sites$Dir,function(x) {if (x=="+") 0 else 1})
    loci.text<-by(sample.sites, 1:nrow(sample.sites), function(x) {if (x$Dir=="+") paste("| ", x$Loc, '>') else paste ('<', x$Loc, " |")})
    loci.alpha<-sapply(sample.sites$SampleAlleleFragments, function(x) {if (x>0) 1 else 0})
    grid.text(
        loci.text,
        unit(sample.sites$Loc,"native"),
        rep(unit(3+sample.sites$Allele,"native"),2),
        hjust=loci.text.just,
        gp=gpar(fontsize=0.5,alpha=loci.alpha)
    )

    # Haplotype sites
    grid.polyline(
        unit(c(sample.sites$MoSite,sample.sites$MoSite),"native"),
        unit(c(rep(3.25,num.alleles),rep(3.75,num.alleles)),"native"),
        id=rep(1:num.alleles,2),
        gp=gpar(lwd=0.5,lineend="round")
    )

    grid.polyline(
        unit(c(sample.sites$FaSite,sample.sites$FaSite),"native"),
        unit(c(rep(2.25,num.alleles),rep(2.75,num.alleles)),"native"),
        id=rep(1:num.alleles,2),
        gp=gpar(lwd=0.5,lineend="round")
    )


    # Fragment lines
    grid.polyline(
        unit(c(sample.sites$MoSite,sample.sites$MoNeighbour),"native"),
        unit(c(rep(3.25,num.alleles),rep(3.75,num.alleles)),"native"),
        id=rep(1:num.alleles,2),
        arrow=arrow(angle=10, length=unit(0.04,"inches"),type="closed"),
        gp=gpar(lwd=0.1,col=rep(sample.sites$Allele,2))
    )

    grid.polyline(
        unit(c(sample.sites$FaSite,sample.sites$FaNeighbour),"native"),
        unit(c(rep(2.25,num.alleles),rep(2.75,num.alleles)),"native"),
        id=rep(1:num.alleles,2),
        arrow=arrow(angle=10, length=unit(0.04,"inches"),type="closed"),
        gp=gpar(lwd=0.1,col=rep(sample.sites$Allele,2))
    )
    
    # Contigs
    
    grid.polyline(
        unit(c(scaffold.agp$Start,scaffold.agp$End),"native"),
        unit(rep(1.5,length(scaffold.agp$Start)*2),"native"),
        id=rep(1:length(scaffold.agp$Start),2),
        arrow=arrow(length=unit(0.02,"inches"),type="closed",ends="both"),
        gp=gpar(lwd=0.5)
    )
    
    # Allele counts
    
    max.allele.count <- max(sample.sites.all$Allele)
    loci.all.dir<-sapply(sample.sites.all$Dir,function(x) {if (x=="+") 1 else -1})

    grid.polyline(
        unit(c(sample.sites.all$Loc-100+200*loci.all.dir,sample.sites.all$Loc+100+200*loci.all.dir),"native"),
        unit(c(0.15+sample.sites.all$Allele/max.allele.count*.8,0.15+sample.sites.all$Allele/max.allele.count*.8),"native"),
        id = rep(1:length(sample.sites.all$Loc),2),
        gp=gpar(lwd=0.01)
    )


    unique.locs<-unique(sample.sites.all$Loc)
    # Max allele count
    grid.polyline(
        unit(c(unique.locs-100,unique.locs+100),"native"),
        unit(rep(0.95,2*length(unique.locs)),"native"),
        id = rep(1:length(unique.locs),2),
        gp=gpar(lwd=0.01)
    )

    # Alleleometer
    grid.polyline(
        unit(c(unique.locs-100,unique.locs-100),"native"),
        unit(c(rep(0.15,length(unique.locs)),rep(0.95,length(unique.locs))),"native"),
        id = rep(1:length(unique.locs),2),
        gp=gpar(lwd=0.01)
    )

    grid.polyline(
        unit(c(unique.locs+100,unique.locs+100),"native"),
        unit(c(rep(0.15,length(unique.locs)),rep(0.95,length(unique.locs))),"native"),
        id = rep(1:length(unique.locs),2),
        gp=gpar(lwd=0.01)
    )
    
    popViewport() # scaffold.vp
    popViewport() # sample.vp
}

plot.scaffold<-function(scaffold.name, scaffold.length, scaffold.agp, scaffold.genotypes, inferred.alleles) {
    grid.newpage()
    # Framing viewport
    pushViewport(viewport(x=0.5,y=0.49,width=0.95,height=0.95, gp=gpar(lineend="butt")))

    # Scaffold master viewport
    scaffold<-viewport(0,0,just=c("left","bottom"),width=1,height=0.1)
    pushViewport(scaffold)

    # Name label viewport
    label<-viewport(0,0,just=c("left","bottom"),width=0.1,height=1)
    pushViewport(label)
    grid.text(scaffold.name,just="left")
    popViewport() #BD.label off

    # scaffold plus axis
    scaffold.line<-viewport(0.1,0,just=c("left","bottom"),width=0.9,height=1,xscale=c(1,scaffold.length))
    pushViewport(scaffold.line)
    grid.lines(unit(c(1,scaffold.length),"native"),c(1,1),gp=gpar(lwd=1,lineend="round"))

    ticks<-seq(0,scaffold.length,by=50000)
    grid.text(sprintf("%7d",ticks), unit(ticks, "native"), 0.7, just="right",rot=90)
    grid.polyline(
    	unit(c(ticks,ticks),"native"),
    	c(rep(0.8,length(ticks)),rep(1,length(ticks))),
    	id=rep(1:length(ticks),2)
    )

    popViewport() # scaffold.line off

    popViewport() # scaffold off

    # samples master viewport
    samples.vp<-viewport(0,0.1,just=c("left","bottom"),width=1,height=0.9)
    pushViewport(samples.vp)

    # Viewports for each individual
    samples<-rev(unique(as.character(scaffold.genotypes$Sample)))
    result<-sapply(samples,function(x) plot.sample(x,inferred.alleles[inferred.alleles$Sample==x & inferred.alleles$Chr==scaffold.name,],samples,scaffold.genotypes,scaffold.length, scaffold.agp))

    popViewport() # samples off
    
}


# Initialise

library(grid)
library(plyr)
library(RColorBrewer)
args<-commandArgs(trailingOnly=T)

scaffold.name<-as.character(args[1])
scaffold.agp<-read.delim(args[2])

inferred.alleles<-read.delim(paste("PstI.all.110802.Hmel1-1_primaryScaffolds.bam.",scaffold.name,".haplotypes.tsv",sep=""))


scaffold.lengths<-read.delim("chr18.scaffolds.tsv",header=F)
names(scaffold.lengths)<-c("Scaffold","Length")
scaffold.length<-scaffold.lengths[scaffold.lengths$Scaffold==scaffold.name,]$Length

names(scaffold.agp)<-c("Scaffold","Start","End","ID","Type","Name","CtgStart","CtgEnd","Orientation","Comment")
scaffold.agp<-scaffold.agp[scaffold.agp$Type == "W" & scaffold.agp$Scaffold==scaffold.name,]

pdf(paste(scaffold.name,".ressites.fragments.pdf",sep=""),width=20,height=25)

genotypes<-read.delim(paste("genotypes.", scaffold.name,".tsv",sep=""))

plot.scaffold(scaffold.name,scaffold.length,scaffold.agp,genotypes,inferred.alleles)

dev.off()

warnings()
