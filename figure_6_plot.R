#!/usr/bin/env Rscript

# figure_6_plot.R
# John W Davey johnomics@gmail.com

library(grid)
library(ggplot2)
library(plyr)
library(scales)
library(RColorBrewer)
args<-commandArgs(trailingOnly=T)

tsv.df.raw<-read.delim(args[1])
reads.df.raw<-read.delim(args[2])
chr18.clean<-read.delim(args[3])

pdf("figure_6_rad_missing_data.pdf",13.30,16)

samples<-c("F1Mo","F1Fa","28","24","73","75")
tsv.df<-tsv.df.raw[tsv.df.raw$Sample %in% samples & tsv.df.raw$FragLen>=400 & tsv.df.raw$FragLen<=10105,]
reads.df<-reads.df.raw[reads.df.raw$Sample %in% samples,]
max.depth<-max(tsv.df$Depth)

tsv.df<-data.frame(tsv.df,
    LogFragLen=log10(tsv.df$FragLen),
    DepthBin=ceiling(tsv.df$Depth/25)/3,
    DepthProp=tsv.df$Depth/max.depth
)

tsv.df<-transform(tsv.df, Sample=factor(tsv.df$Sample,levels=samples,labels=c("Mother","Father","AB","AC","CB","CC")))
reads.df<-transform(reads.df, Sample=factor(reads.df$Sample,levels=samples,labels=c("Mother","Father","AB","AC","CB","CC")))

chr18.clean.select<-chr18.clean[chr18.clean$Sample %in% samples & chr18.clean$AlleleCount>0 & (!is.na(chr18.clean$MoFragLen) | !is.na(chr18.clean$FaFragLen) | !is.na(chr18.clean$SharedFragLen)),]
chr18.clean.select<-transform(chr18.clean.select,Sample=factor(chr18.clean.select$Sample,levels=c("F1Mo","F1Fa","28","24","73","75"), labels=c("Mother (AC)", "Father (BC)", "AB","AC","CB","CC")))


allele.col<-data.frame(Allele=c("1","2","3"),Colour=c("black","red","blue"))

# Master viewport
pushViewport(viewport(x=0.5,y=0.49,width=0.95,height=0.95,gp=gpar(fontfamily="Palatino",lineend="butt")))

alleles.vp<-viewport(0,0,just=c("left","bottom"),width=0.88,height=0.8)
pushViewport(alleles.vp)



xscale.vp<-viewport(0.075,0,just=c("left","bottom"),width=0.925,height=0.1,xscale=c(1,length(reads.df$Sample)*10+0.1))
pushViewport(xscale.vp)

grid.lines(
    x=c(1,length(reads.df$Sample)*10),
    y=1,
    default.units="native"
)


# Read depth scales
grid.polyline(
    x=c(seq(1.9,length(reads.df$Sample)*10,10),seq(5.9,length(reads.df$Sample)*10,10)),
    y=rep(0.9,2*length(reads.df$Sample)),
    id=rep(1:length(reads.df$Sample),2),
    default.units="native",
    gp=gpar(lineend="round")
)

grid.polyline(
    x=c(seq(6.1,length(reads.df$Sample)*10+0.1,10),seq(10.1,length(reads.df$Sample)*10+0.1,10)),
    y=rep(0.9,2*length(reads.df$Sample)),
    id=rep(1:length(reads.df$Sample),2),
    default.units="native",
    gp=gpar(lineend="round")
)


# Read depth tick lines
grid.polyline(
    x=c(seq(1.9,length(reads.df$Sample)*10,10),seq(1.9,length(reads.df$Sample)*10,10)),
    y=c(rep(0.9,length(reads.df$Sample)),rep(0.85,length(reads.df$Sample))),
    id=rep(1:length(reads.df$Sample),2),
    default.units="native"
)

grid.polyline(
    x=c(seq(3.9,length(reads.df$Sample)*10,10),seq(3.9,length(reads.df$Sample)*10,10)),
    y=c(rep(0.9,length(reads.df$Sample)),rep(0.85,length(reads.df$Sample))),
    id=rep(1:length(reads.df$Sample),2),
    default.units="native"
)

grid.polyline(
    x=c(seq(5.9,length(reads.df$Sample)*10,10),seq(5.9,length(reads.df$Sample)*10,10)),
    y=c(rep(0.9,length(reads.df$Sample)),rep(0.85,length(reads.df$Sample))),
    id=rep(1:length(reads.df$Sample),2),
    default.units="native"
)

grid.polyline(
    x=c(seq(6.1,length(reads.df$Sample)*10+0.1,10),seq(6.1,length(reads.df$Sample)*10+0.1,10)),
    y=c(rep(0.9,length(reads.df$Sample)),rep(0.85,length(reads.df$Sample))),
    id=rep(1:length(reads.df$Sample),2),
    default.units="native"
)

grid.polyline(
    x=c(seq(8.1,length(reads.df$Sample)*10+0.1,10),seq(8.1,length(reads.df$Sample)*10+0.1,10)),
    y=c(rep(0.9,length(reads.df$Sample)),rep(0.85,length(reads.df$Sample))),
    id=rep(1:length(reads.df$Sample),2),
    default.units="native"
)

grid.polyline(
    x=c(seq(10.1,length(reads.df$Sample)*10+0.1,10),seq(10.1,length(reads.df$Sample)*10+0.1,10)),
    y=c(rep(0.9,length(reads.df$Sample)),rep(0.85,length(reads.df$Sample))),
    id=rep(1:length(reads.df$Sample),2),
    default.units="native"
)

grid.text(
    70,
    x=c(seq(1.9,length(reads.df$Sample)*10,10)),
    y=0.80,
    just=c("centre","top"),
    default.units="native"
)

grid.text(
    35,
    x=c(seq(3.9,length(reads.df$Sample)*10,10)),
    y=0.80,
    just=c("centre","top"),
    default.units="native"
)

grid.text(
    0,
    x=c(seq(6,length(reads.df$Sample)*10,10)),
    y=0.80,
    just=c("centre","top"),
    default.units="native"
)

grid.text(
    35,
    x=c(seq(8.1,length(reads.df$Sample)*10,10)),
    y=0.80,
    just=c("centre","top"),
    default.units="native"
)

grid.text(
    70,
    x=c(seq(10.1,length(reads.df$Sample)*10+0.1,10)),
    y=0.80,
    just=c("centre","top"),
    default.units="native"
)

grid.text(c("Mother","Father","AB","AC","CB","CC"),seq(6,length(reads.df$Sample)*10,10),0.1,default.units="native")
grid.text(c("A","B","A","A","C","C"),seq(4,length(reads.df$Sample)*10,10),0.3,default.units="native")
grid.text(c("C","C","B","C","B","C"),seq(8,length(reads.df$Sample)*10,10),0.3,default.units="native")
grid.text("Reads",seq(6,length(reads.df$Sample)*10,10),0.5,default.units="native")

popViewport() #xscale.vp

yscale.vp<-viewport(0,0.1,just=c("left","bottom"),width=0.075,height=0.9,yscale=c(log10(400),log10(10500)))
pushViewport(yscale.vp)

grid.text(
    "Restriction Fragment Length (kbp)",
    x =0.1,
    rot=90
)

grid.text(
    c("0.4","1","2","3","4","5","6","7","8","9","10"),
    x=0.75,
    just="right",
    y=log10(c(400,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)),
    default.units="native"
)
grid.polyline(
    x=rep(c(0.8,1),11),
    y=rep(log10(c(400,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)),2),
    id=rep(1:11,2),
    default.units="native"
)

grid.lines(
    x=1,
    y=log10(c(400,10000)),
    default.units="native"
)

popViewport() #yscale.vp

samples.vp<-viewport(0.075,0.1,just=c("left","bottom"),width=0.925,height=0.9,xscale=c(1,length(reads.df$Sample)*10),yscale=c(log10(400),log10(10500)))
pushViewport(samples.vp)


     vjust<-sapply(tsv.df$Dir,    function(x) {if (x=="+")  1 else -1})
sample.pos<-sapply(tsv.df$Sample, function(s) {which(reads.df$Sample==s)})*10-10

grid.polyline(
    x = c(sample.pos+6+tsv.df$HapDir*0.1,sample.pos+6+tsv.df$HapDir*0.1+tsv.df$HapDir*tsv.df$DepthProp*4),
    y = c(tsv.df$LogFragLen+0.001*vjust,tsv.df$LogFragLen+0.001*vjust),
    id = rep(1:length(tsv.df$Depth),2),
    default.units="native",
    gp = gpar(lwd=0.25,col=rep(sapply(tsv.df$HapCol,function(x){as.character(allele.col[allele.col$Allele==x,]$Colour)}),2))
)
popViewport() # samples.vp
popViewport() # alleles.vp

missing.vp<-viewport(0.9,0,just=c("left","bottom"),width=0.1,height=0.8)
pushViewport(missing.vp)

xscale.vp<-viewport(0,0,just=c("left","bottom"),width=1,height=0.1,xscale=c(0,length(reads.df.raw$Sample)))
pushViewport(xscale.vp)
grid.lines(
    x=c(0,length(reads.df.raw$Sample)),
    y=1,
    default.units="native"
)
ind.ticks<-seq(0,length(reads.df.raw$Sample),5)

grid.polyline(
    x=rep(ind.ticks,2),
    y=c(rep(0.9,length(ind.ticks)),rep(1,length(ind.ticks))),
    id=rep(1:length(ind.ticks),2),
    default.units="native"
)
ind.tick.labels<-seq(0,length(reads.df.raw$Sample),15)
grid.text(
    ind.tick.labels,
    x=ind.tick.labels,
    y=0.8,
    just=c("centre","top"),
    default.units="native"
)

grid.text("Individuals",0.5,0.5,default.units="npc")

popViewport() #xscale.vp

individuals.vp<-viewport(0,0.1,just=c("left","bottom"),width=1,height=0.9,xscale=c(0,length(reads.df.raw$Sample)),yscale=c(log10(400),log10(10500)))
pushViewport(individuals.vp)

grid.polyline(
    x = c(rep(0,length(tsv.df$Missing)),tsv.df$Missing),
    y = c(tsv.df$LogFragLen+0.0025*vjust,tsv.df$LogFragLen+0.0025*vjust),
    id = rep(1:length(tsv.df$Missing),2),
    default.units="native",
    gp=gpar(lwd=0.25,col="black")
)
grid.polyline(
    x = c(tsv.df$Missing,tsv.df$Missing+tsv.df$HetHom),
    y = c(tsv.df$LogFragLen+0.0025*vjust,tsv.df$LogFragLen+0.0025*vjust),
    id = rep(1:length(tsv.df$Missing),2),
    default.units="native",
    gp=gpar(lwd=0.25,col="red")
)

grid.lines(x=0, y=c(log10(400),log10(10500)),default.units="native",gp=gpar(lwd=1,col="gray"))
grid.lines(x=length(reads.df.raw$Sample), y=c(log10(400),log10(10500)),default.units="native",gp=gpar(lwd=1,col="gray"))

popViewport() # individuals.vp
popViewport() # missing.vp

allelecount.vp<-viewport(0.015,0.8,just=c("left","bottom"),width=0.96,height=0.2)
print(
    ggplot(chr18.clean.select,
           aes(SharedFragLen,SampleAlleleReads,colour=factor(AlleleCount)))+
           facet_grid(.~Sample)+
           geom_point(size=1)+
           ylim(c(0,100))+
           xlim(c(1,1200000))+
           scale_x_continuous(trans=log10_trans(),breaks=c(400,1000,4000,10000,40000),labels=c("0.4","1","4","10","40"))+
           theme_bw(base_family="Palatino")+
           scale_colour_manual(name="Allele\nCopies",values=c("mediumorchid1","forestgreen"))+
           xlab("Restriction Fragment Length (kbp)")+
           ylab("Read Depth"),
    vp=allelecount.vp
)

popViewport() # Master
dev.off()
warnings()
