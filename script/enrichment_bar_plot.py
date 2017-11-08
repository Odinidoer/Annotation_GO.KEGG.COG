#!/usr/bin/env python
#coding:utf-8
#fengyitong
#20171010

import argparse
import os
import sys
sys.path.append('/mnt/ilustre/users/hui.wan/.local/lib/python2.7/site-packages/Mako-1.0.6-py2.7.egg') 
from mako.template import Template

parser = argparse.ArgumentParser(
    description="""
Usage:        perl enrich_barplot.pl [options]
Description:  This program is used for ploting GO enrichment or KEGG enrichment analysis results.
Contact:      yitong.feng@majorbio.com  Time: 2017.10
Version:      1.0
""")
parser.add_argument(
    "-i",
    dest="file",
	required=True,
    type=str,
    help="please input file")
parser.add_argument(
    "-w",
    dest="w",
    type=str,
    default = '16',
    help='plot width ,defalt: 16')
parser.add_argument(
    "-ph",
    dest="ph",
    type=str,
    default = '7',
    help='plot heigth ,defalt: 7')
parser.add_argument(
    "-pc",
    dest="pc",
    type=str,
    default = '0.9',
    help='cex of the pathway names in the x axis,default: 0.9')
parser.add_argument(
    "-sc",
    dest="sc",
    type=str,
    default = '0.3',
    help='cex of the stars in the top of each bars,default:0.3')
parser.add_argument(
    "-mb",
    dest="mb",
    type=str,
    default = '13',
    help='white space of bottom,default:13')
parser.add_argument(
    "-ml",
    dest="ml",
    type=str,
    default = '6',
    help='white space of left,default:6')
parser.add_argument(
    "-o",
    dest="o",
    type=str,
	required=True,
    help='out_name of enrich_barplot')	
args = parser.parse_args()
Rcmd='''options(warn=-1)
w<-${w}
h<-${ph}
mb<-${mb}
ml<-${ml}
## define colors
getMyColor<-function(col1="#008B45",col2="white",cutn=10010,targets){
	targets=round(targets,4)
	ramp_cc<-colorRamp(as.vector(c(col1,col2)))
	continue_cols<-paste(rgb(ramp_cc(round(seq(0,1,length=cutn),4)),max=255),"E5",sep="")
	num2col<-as.data.frame(cbind(round(seq(0,1,length=cutn),4),continue_cols))
	names(num2col)<-c("number","color")
	targetcolor<-vector()
	for(i in 1:length(targets)){
		targetcolor[i]<-as.character(num2col[which(num2col[[1]]==targets[i]),2])
	}
	return(targetcolor)
}
PATHWAY_enrich_type<-read.delim("${i}",header=T,sep="\t",check.names=F)
if(nrow(PATHWAY_enrich_type) < 15){
	print("The number of enriched KEGG Pathways is less than 15,the picture would not been shown!")
	}else{
		### color card:
		ramp_cc<-colorRamp(as.vector(c("#008B45","white")))
		continue_cols<-paste(rgb(ramp_cc(round(seq(0,1,length=10010),4)),max=255),"E5",sep="")
		PATHWAY_enrich_raw<-PATHWAY_enrich_type[order(PATHWAY_enrich_type[,6]),]
        ### top bars to be view
		if(nrow(PATHWAY_enrich_raw)>60){             
			PATHWAY_enrich<-PATHWAY_enrich_raw[1:60,]
		}else{PATHWAY_enrich<-PATHWAY_enrich_raw}
		RichFactor1<-sub("\\\|",";",PATHWAY_enrich[,5])
		PATHWAY_enrich<-PATHWAY_enrich[order(as.numeric(PATHWAY_enrich[,3]),as.numeric(PATHWAY_enrich[,2])),]
        RichFactor<-NA
        for (i in RichFactor1){RichFactor<-c(RichFactor,as.numeric(unlist(strsplit(i,split = ";"))[1])/as.numeric(unlist(strsplit(i,split = ";"))[2]))}	
        RichFactor<-round(RichFactor[-1],2)		
		### significant stars:
		Pvalue<-as.numeric(PATHWAY_enrich[[6]]) ## pvalue
		n_3tars<-which(Pvalue<=0.001)
		n_2tars<-setdiff(which(Pvalue<=0.01),n_3tars)
		n_1tars<-setdiff(which(Pvalue<=0.05),c(n_3tars,n_2tars))
		barcols<-getMyColor(targets=Pvalue)
		### pathway class:
		name2typeI<-PATHWAY_enrich[,c(2,3)]
		shortName<-gsub("Organismal Systems","OS",gsub("Human Diseases","HD",gsub("Genetic Information Processing","GIP",gsub("Environmental Information Processing","EIP",gsub("Drug Development","DD",gsub("Cellular Processes","CP",gsub("Metabolism","M",gsub("biological_process","BP",gsub("cellular_component","CC",gsub("molecular_function","MF",name2typeI[[2]]))))))))))
		xlabs<-paste(PATHWAY_enrich[[2]],shortName,sep=":  ")
		
		### ploting
		pdf('${o}',width=w,height=h)
		skip<-max(RichFactor)/50
		par(fig=c(0,0.92,0,1),mar = c(${mb},${ml},4,0.1)+0.1,new=F)
		bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
		text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=${pc},font=1)
		
		#text(bar,rep(-0.01,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=${pc},font=1)
		if(length(n_3tars)!=0){text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=${sc},col="black")}
		if(length(n_2tars)!=0){text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=${sc},col="black")}
		if(length(n_1tars)!=0){text(bar[n_1tars],RichFactor[n_1tars]+skip,"*",cex=${sc},col="black")}
		
		### legend for bar colors:
		par(fig=c(0.91,1,0,0.9),mar = c(${mb},0.1,4,4)+0.1,new=T)
		Bars<-barplot(rep(0.5,10000),horiz=T,col=continue_cols[10000:1],border=continue_cols[10000:1],space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
		text(rep(0.6,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
		text(0,length(continue_cols)+500,"Pvalue",adj=0,xpd=T)
		dev.off()}'''

f = Template(Rcmd)
Rcmd = f.render(i=args.file,w=args.w,ph=args.ph,pc=args.pc,sc=args.sc,mb=args.mb,ml=args.ml,o=args.o)
fout=open('%s.cmd.r' %(args.o),'w')	
fout.write(Rcmd)
fout.close()
cmd = '/mnt/ilustre/users/shuirong.zhang/softwares/bin/Rscript %s.cmd.r' %(args.o)
os.system(cmd)
