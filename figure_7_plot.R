#!/usr/bin/env Rscript

# figure_7_plot.R
# John W Davey johnomics@gmail.com

library(ggplot2)
library(plyr)
library(RColorBrewer)
args<-commandArgs(trailingOnly=T)

tool.df<-read.delim(args[1])
tool.df.reads<-transform(tool.df, Tool=factor(tool.df$Tool,levels=c("Reads")))[tool.df$Tool == "Reads",]
tool.df.tools<-transform(tool.df, Tool=factor(tool.df$Tool,levels=c("RADtools Clusters","GATK","Stacks Clusters","Stacks Genotypes")))[tool.df$Tool != "Reads",]

readprop.pdf<-args[2]
fragskew.pdf<-args[3]
if(is.na(readprop.pdf) | is.na(fragskew.pdf)) stop("Please specify two filenames for two PDFs: one for read proportions, one for tool performance")

breaks<-c(0,0.5,1)

pdf(readprop.pdf,7,6.2)
ggplot(tool.df.reads,
    aes(FragLen1,FragLen2,colour=Proportion))+
    geom_point(size=3)+
    scale_x_log10(limits=c(500,60000),breaks=c(500,1000,2500,5000,10000,25000,50000))+
    scale_y_log10(limits=c(500,60000),breaks=c(500,1000,2500,5000,10000,25000,50000))+
    scale_colour_gradientn(
        name="Read\nProportion\n",
        colours=rev(brewer.pal(3,"Set1")),
        breaks=c(0,0.5,1),
        labels=format(breaks),
        limits=c(0,1)
    )+
    xlab("Restriction Fragment Length 1")+
    ylab("Restriction Fragment Length 2")+
    facet_wrap(~Tool)+
    theme_bw(base_family="Palatino")

dev.off()

pdf(fragskew.pdf,7,6.2)
ggplot(tool.df.tools,
    aes(FragLen1,FragLen2,colour=Proportion))+
    geom_point()+
    scale_x_log10(limits=c(500,60000),breaks=c(500,1000,2500,5000,10000,25000,50000))+
    scale_y_log10(limits=c(500,60000),breaks=c(500,1000,2500,5000,10000,25000,50000))+
    scale_colour_gradientn(
        name="Error\nProportion\n",
        colours=rev(brewer.pal(3,"Set1")),
        breaks=breaks,
        labels=format(breaks),
        limits=c(0,1)
    )+
    xlab("Restriction Fragment Length 1")+
    ylab("Restriction Fragment Length 2")+
    facet_wrap(~Tool,ncol=2)+
    theme_bw(base_family="Palatino")+
    opts(axis.text.x = theme_text(size = 8,family="Palatino"),axis.text.y = theme_text(size = 8,family="Palatino"))
    
dev.off()

warnings()
