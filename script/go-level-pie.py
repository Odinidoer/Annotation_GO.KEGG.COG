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
go_level2 to pie picture
""")

parser.add_argument(
    "-i",
    dest="infile",
	required=True,
    type=str,
    help="please input file")
parser.add_argument(
    "-o",
    dest="outfile",
	required=True,
    type=str,
    help="the name of  outfile")
parser.add_argument(
    "-lv",
    dest="lv",
    type=str,
    default = '2',
    help='the number of go_level ,defalt: 2')
parser.add_argument(
    "-w",
    dest="w",
    type=str,
    default = '10',
    help='plot width ,defalt: 10')
parser.add_argument(
    "-ph",
    dest="ph",
    type=str,
    default = '8',
    help='plot heigth ,defalt: 8')
parser.add_argument(
    "-pc",
    dest="pc",
    type=str,
    default = '0.9',
    help='cex of the  labels of pies in the x axis,default: 0.9')
parser.add_argument(
    "-lc",
    dest="lc",
    type=str,
    default = '1',
    help='cex of the labels of go_level2 in the topleft,default:1')
	
args = parser.parse_args()

Rcmd='''library(RColorBrewer)
lv<-${lv}
w<-${w}
h<-${ph}
mb<-${pc}
ml<-${lc}
GO <- read.delim("${i}",sep = "\t",header = TRUE)
GO2 <- GO[GO$level == lv,]
x<-c(0.01,0.31,0.32,0.62,0.63,0.93)
xl<-c(0.08,0.38,0.40,0.68,0.69,1)
y<-c(0.2,0.5,0.5,0.73)
ramp_bp<-colorRamp(as.vector(c("#6B8E23","white")))
heatcol_BP2<-paste(rgb(ramp_bp(seq(0,1,length=10)),max=255),"E5",sep="")
ramp_cc<-colorRamp(as.vector(c("#7EC0EE","white")))
heatcol_CC2<-paste(rgb(ramp_cc(seq(0,1,length=10)),max=255),"E5",sep="")
ramp_mf<-colorRamp(as.vector(c("#FF69B4","white")))
heatcol_MF2<-paste(rgb(ramp_mf(seq(0,1,length=10)),max=255),"E5",sep="")    
ColorList<-brewer.pal(10,"Paired")

pdf('${o}',w,h)
layout(matrix(1:6,2,3,byrow=F),heights=c(3,2))

BP_2<-GO2[which(GO2$namespace=="biological_process"),]
BP_2_sort<-BP_2[order(BP_2$num_of_Accs,decreasing=T),][1:ifelse(nrow(BP_2)>10,10,nrow(BP_2)),]
BP_2_splice<-as.numeric(as.character(BP_2_sort$num_of_Accs))
BP_2_label<-as.character(BP_2_sort$num_of_Accs)
BP_2_legend<-paste(as.character(BP_2_sort$name)," (",BP_2_label,")",sep="")

CC_2<-GO2[which(GO2$namespace=="cellular_component"),]
CC_2_sort<-CC_2[order(CC_2$num_of_Accs,decreasing=T),][1:ifelse(nrow(CC_2)>10,10,nrow(CC_2)),]
CC_2_splice<-as.numeric(as.character(CC_2_sort$num_of_Accs))
CC_2_label<-as.character(CC_2_sort$num_of_Accs)
CC_2_legend<-paste(as.character(CC_2_sort$name)," (",CC_2_label,")",sep="")

MF_2<-GO2[which(GO2$namespace=="molecular_function"),]
MF_2_sort<-MF_2[order(MF_2$num_of_Accs,decreasing=T),][1:ifelse(nrow(MF_2)>10,10,nrow(MF_2)),]
MF_2_splice<-as.numeric(as.character(MF_2_sort$num_of_Accs))
MF_2_label<-as.character(MF_2_sort$num_of_Accs)
MF_2_legend<-paste(as.character(MF_2_sort$name)," (",MF_2_label,")",sep="")

par(mai=c(0,0,0,0))
par(fig=c(x[1],x[2],y[3],y[4]),new=F)
pie(MF_2_splice,labels=MF_2_label,col=heatcol_MF2,clockwise=T,cex=mb,radius=0.89,border="white")
legend("topright","MF",bty="n",cex=2,text.col="#FF69B4",text.font=2)
par(fig=c(xl[1],xl[2],y[1],y[2]),new=T)
plot(1:2,xlab="",ylab="",axes=F,type="n")
legend("topleft",legend=MF_2_legend,fill=heatcol_MF2,cex=ml,bty="n")

par(fig=c(x[3],x[4],y[3],y[4]),new=T)
pie(CC_2_splice,labels=CC_2_label,col=heatcol_CC2,clockwise=T,cex=mb,radius=0.89,border="white")
legend("topright","CC",bty="n",cex=2,text.col="#7EC0EE",text.font=2)
par(fig=c(xl[3],xl[4],y[1],y[2]),new=T)
plot(1:2,xlab="",ylab="",axes=F,type="n")
legend("topleft",legend=CC_2_legend,fill=heatcol_CC2,cex=ml,bty="n")

par(fig=c(x[5],x[6],y[3],y[4]),new=T)
pie(BP_2_splice,labels=BP_2_label,col=heatcol_BP2,clockwise=T,cex=mb,radius=0.89,border="white")
legend("topright","BP",bty="n",cex=2,text.col="#6B8E23",text.font=2)
par(fig=c(xl[5],xl[6],y[1],y[2]),new=T)
plot(1:2,xlab="",ylab="",axes=F,type="n")
legend("topleft",legend=BP_2_legend,fill=heatcol_BP2,cex=ml,bty="n")'''

f = Template(Rcmd)
Rcmd = f.render(i=args.infile,o=args.outfile,lv=args.lv,w=args.w,ph=args.ph,pc=args.pc,lc=args.lc)
fout=open('%s.r' %args.infile,'w')	
fout.write(Rcmd)
fout.close()
os.system('Rscript %s.r' %args.infile)