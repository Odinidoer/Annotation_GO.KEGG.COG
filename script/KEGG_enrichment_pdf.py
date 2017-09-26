#!/usr/bin/env python
#coding:utf-8
#yanjun
#20170918

import argparse
import ConfigParser
try:
	from mako.template import Template
except:
	import sys
	sys.path.append('/mnt/ilustre/users/hui.wan/.local/lib/python2.7/site-packages/Mako-1.0.6-py2.7.egg')
	from mako.template import Template
else:
	exit('cannot import mako')

parser=argparse.ArgumentParser(description="plot KEGG enrichment pdf")
parser.add_argument("-config",type=str,required=True,help="config.ini")
parser.add_argument("-enrichmentxls",type=str,required=True,help="KEGG_enrichment.xls for plot")
parser.add_argument("-w",type=int,help="plot width ,design by num_of_accs automaticly:16")
parser.add_argument("-h",type=int,help="plot width ,design by length_of_accs`s_name automaticly:7")
parser.add_argument("-pc",type=float,default=0.9,help=" cex of the pathway names in the x axis,default: 0.9")
parser.add_argument("-lc",type=float,default=1,help="cex of the class names in the topleft,default:1")
parser.add_argument("-sc",type=float,default=0.3,help="cex of the stars in the top of each bars,default:0.3")
parser.add_argument("-mb",type=float,default=13,help="white space of bottom,default:13")
parser.add_argument("-ml",type=float,default=6,help="white space of left,default:6")
parser.add_argument("-out_plot",type=str,required=True,help="outfile for plot")

args=parser.parse_args()

cmd_r = '''
options(warn=-1)
#w<-$opts{w}
#h<-$opts{h}
mb<-$mb
ml<-$ml

## define colors
getMyColor<-function(col1="#008B45",col2="white",cutn=10010,targets){
	targets=round(targets,4) ## round 保留小数点位数
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

PATHWAY_type<-read.delim(${index_file}",sep="\t",header=T,check.names=F)
PATHWAY_enrichraw<-read.delim(${enrichmentxls_file},header=T,sep="\t",check.names=F)
PATHWAY_enrich_type<-merge(PATHWAY_enrichraw,PATHWAY_type,by="#Term")[,c("#Term","ratio","pvalue","typeII","typeI")]

### color card:
ramp_cc<-colorRamp(as.vector(c("#008B45","white")))
continue_cols<-paste(rgb(ramp_cc(round(seq(0,1,length=10010),4)),max=255),"E5",sep="")
PATHWAY_enrich_raw<-PATHWAY_enrich_type[order(PATHWAY_enrich_type[,5],as.numeric(PATHWAY_enrich_type[,4]),
as.numeric(PATHWAY_enrich_type[,3]),as.numeric(PATHWAY_enrich_type[,2])),]
### top bars to be view
if(nrow(PATHWAY_enrich_raw)>60){
	topNs<-which(PATHWAY_enrich_raw[[3]]<0.5)## pvalue >0.5 were ignored 
	PATHWAY_enrich<-PATHWAY_enrich_raw[topNs,]
	}
else{
	PATHWAY_enrich<-PATHWAY_enrich_raw
	}

RichFactor<-round(PATHWAY_enrich[,2],2)            
### significant stars:
Pvalue<-as.numeric(PATHWAY_enrich[[3]]) ## pvalue
n_3tars<-which(Pvalue<=0.001)
n_2tars<-setdiff(which(Pvalue<=0.01),n_3tars)
n_1tars<-setdiff(which(Pvalue<=0.05),c(n_3tars,n_2tars))
barcols<-getMyColor(targets=Pvalue)
### pathway class:
name2typeI<-PATHWAY_enrich[,c(1,6)]
shortName<-gsub("Organismal Systems","OS",gsub("Human Diseases","HD",gsub("Genetic Information Processing",
"GIP",gsub("Environmental Information Processing","EIP",gsub("Drug Developmen","DD",gsub("Cellular Processes","CP",gsub("Metabolism","M",name2typeI[[2]])))))))
xlabs<-paste(PATHWAY_enrich[[1]],shortName,sep=":   ")
typeII<-c("KEGG Pathway Class:","EIP: Environmental Information Processing","GIP: Genetic Information Processing",
"CP : Cellular Processes","OS : Organismal Systems","DD : Drug Development","HD : Human Diseases","M  : Metabolism")
### ploting
pdfname <-"${pdf_name}.pathway.pdf"
pdf(file=pdfname,width=w,height=h)
skip<-max(RichFactor)/50

par(fig=c(0,0.92,0,1),mar = c($opts{mb},$opts{ml},4,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,
ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=$opts{pc},font=1)

if(length(n_3tars)!=0){text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=$opts{sc},col="black")}
if(length(n_2tars)!=0){text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=$opts{sc},col="black")}
if(length(n_1tars)!=0){text(bar[n_1tars],RichFactor[n_1tars]+skip,"*",cex=$opts{sc},col="black")}
legend("topleft",typeII,cex=$opts{lc},bty="n",inset=$opts{is})

### legend for bar colors:
par(fig=c(0.91,1,0,0.9),mar = c($opts{mb},0.1,4,4)+0.1,new=T)
Bars<-barplot(rep(0.5,10000),horiz=T,col=continue_cols[10000:1],border=continue_cols[10000:1],
space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.6,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+500,"Pvalue",adj=0,xpd=T)
dev.off()''