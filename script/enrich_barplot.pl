#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions (\%opts,"t=s","i=s","w=i","h=i","pc=f","lc=f","sc=f","mb=f","ml=f","is=f");
##
my $usage = <<"USAGE";
Usage:        perl enrich_barplot.pl [options]
Description:  This program is used for ploting GO enrichment or KEGG enrichment analysis results.
              enrichment file name format:
	      KEGG: *kegg_enrichment.xls
	      GO  : *enrichment.detail.xls
Contact:      yan.wang\@majorbio.com  Time: 2014.06.25
Version:      1.0
Options:
       input file related:
			-i		FILE     GO or KEGG enrichment result table[required]
					(must be obtained from find_enrichment.py or exact_goatools.pl)
			-t		STRING   enrichment type: GO-PATHWAY-DISEASE
       plot related:
			-w		INT      plot width ,defalt: 16
			-h		INT      plot width ,defalt: 7
			-pc		FLOAT    cex of the pathway names in the x axis,default: 0.9
			-lc		FLOAT    cex of the class names in the topleft,default:1
			-sc		FLOAT    cex of the stars in the top of each bars,default:0.3
			-mb		FLOAT	white space of bottom,default:13
			-ml		FLOAT	white space of left,default:6
			-is		FlLOAT	inset percent of the legend,default 0
USAGE

die $usage if ( !defined $opts{i});
die $usage if ( !defined $opts{t});

#define defaults
$opts{w}=defined$opts{w}?$opts{w}:16;
$opts{h}=defined$opts{h}?$opts{h}:7;
$opts{pc}=defined$opts{pc}?$opts{pc}:0.9;
$opts{lc}=defined$opts{lc}?$opts{lc}:1;
$opts{sc}=defined$opts{sc}?$opts{sc}:0.8;
$opts{mb}=defined$opts{mb}?$opts{mb}:13;
$opts{ml}=defined$opts{ml}?$opts{ml}:6;
$opts{is}=defined$opts{is}?$opts{is}:0;


if(lc($opts{t}) eq "go")
{
	open GO, $opts{i} or die $!;
	open GOO, ">$opts{i}.go" or die $!;
	while(<GO>)
	{
		chomp;
		my @tmp = split /\t/;
		if(/^id/)
		{
			print GOO "description\tratio\tp_bonferroni\ttype\tpvalue\n";
		}else{
			my $fenzi = (split /\//, $tmp[3])[0];
			my $fenmu = (split /\//, $tmp[4])[0];
			my $ratio = $fenzi / $fenmu;
			next if($tmp[1] =~  /p/);
#			my @num = split / /, $tmp[2];
#			my $name = $tmp[2];
#			if(@num > 5)
#			{
#				$name = join " ", @num[0..4], "\b...";
#			}else{
#				
#			}
			print GOO "$tmp[2]\t$ratio\t$tmp[6]\t$tmp[7]\t$tmp[5]\n";
		}
	}
	close GO;
	close GOO;
	system("head -81 $opts{i}.go > $opts{i}.go.draw");

}

if(lc($opts{t}) eq "pathway" || lc($opts{t}) eq "disease"){
	open EN, $opts{i} or die $!;
	open PWP, ">$opts{i}.pathway" or die $!;
	open PWD, ">$opts{i}.disease" or die $!;
	open PWP_1, ">$opts{i}.pathway.xls" or die $!;
	open PWD_1, ">$opts{i}.disease.xls" or die $!;
	
	$/ = "\n#";
	<EN>;<EN>;<EN>;
	my $p = <EN>;
	my @lines = split /\n/, $p;
	print PWP "#Term\tratio\tCorrected P-Value\tpvalue\n";
	print PWP_1 "#Term\tDatabase\tId\tSample number\tBackground number\tP-Value\tCorrected P-Value\tProteins\tHyperlink\n";	
	for(my $i = 1; $i < @lines; $i ++)
	{
		next if($lines[$i] !~ /^\w+/);
		my @tmp = split /\t/, $lines[$i];
		my $ratio = $tmp[3] / $tmp[4];
#		my @num = split / /, $tmp[0];
#		my $name = $tmp[0];
#		if(@num > 5)
#		{
#			$name = join " ", @num[0..4], "\b...";
#		}else{
#			
#		}
		print PWP "$tmp[0]\t$ratio\t$tmp[6]\t$tmp[5]\n";
		print PWP_1 $lines[$i]."\n";
		
	}
	my $d = <EN>;
	my @line = split /\n/, $d;
	print PWD "#Term\tratio\tCorrected P-Value\tpvalue\n";
	print PWD_1 "#Term\tDatabase\tId\tSample number\tBackground number\tP-Value\tCorrected P-Value\tProteins\tHyperlink\n";
	for(my $i = 1; $i < @line; $i ++)
	{
		next if($line[$i] !~ /^\w+/);
		my @tmp = split /\t/, $line[$i];
		my $ratio = $tmp[3] / $tmp[4];
#		my @num = split / /, $tmp[0];
#		my $name = $tmp[0];
#		if(@num > 5)
#		{
#			$name = join " ", @num[0..4], "\b...";
###		}
		print PWD "$tmp[0]\t$ratio\t$tmp[6]\t$tmp[5]\n";
		print PWD_1 $lines[$i]."\n";
	}
	close EN;
	close PWP;
	close PWD;
	close PWP_1;
	close PWD_1;
	## add pathway class
	system "perl /mnt/ilustre/users/ting.kuang/scripts/proteomics/pipeline/bin/tabletools_add.pl -i /mnt/lustre-2/users/wangyan/rna/database/KEGG/pathwayName_typeII_typeI -t $opts{i}.pathway.xls -n 1 -headt T -headi T > $opts{i}.pathway_addclass.xls";
	system "perl /mnt/ilustre/users/ting.kuang/scripts/proteomics/pipeline/bin/tabletools_add.pl -i /mnt/lustre-2/users/wangyan/rna/database/KEGG/pathwayName_typeII_typeI -t $opts{i}.disease.xls -n 1 -headt T -headi T > $opts{i}.disease_addclass.xls";
}
if(!(lc($opts{t}) eq "go" || lc($opts{t}) eq "pathway" || lc($opts{t}) eq "disease")){
	die "please choose correct type: GO or PATHWAY or DISEASE, case-insensitive!!!";
}

open RCMD, ">$opts{i}.r";
print RCMD "
options(warn=-1)
w<-$opts{w}
h<-$opts{h}
mb<-$opts{mb}
ml<-$opts{ml}
type<-\"$opts{t}\"
mycol <-c(119,132,147,454,89,404,123,529,463,104,552,28,54,84,256,100,558,43,652,31,610,477,588,99,81)
## define colors
getMyColor<-function(col1=\"#008B45\",col2=\"white\",cutn=10010,targets){
	targets=round(targets,4)
	ramp_cc<-colorRamp(as.vector(c(col1,col2)))
	continue_cols<-paste(rgb(ramp_cc(round(seq(0,1,length=cutn),4)),max=255),\"E5\",sep=\"\")
	num2col<-as.data.frame(cbind(round(seq(0,1,length=cutn),4),continue_cols))
	names(num2col)<-c(\"number\",\"color\")
	targetcolor<-vector()
	for(i in 1:length(targets)){
		targetcolor[i]<-as.character(num2col[which(num2col[[1]]==targets[i]),2])
	}
	return(targetcolor)
}

if(type==\"PATHWAY\"){
	PATHWAY_type<-read.delim(\"/mnt/ilustre/users/ting.kuang/scripts/proteomics/pipeline/bin/kegg/pathwayId_pathwayName_typeII_typeI\",sep=\"\\t\",header=T,check.names=F)
    PATHWAY_enrichraw<-read.delim(\"$opts{i}.pathway\",header=T,sep=\"\\t\",check.names=F)
	PATHWAY_enrich_type<-merge(PATHWAY_enrichraw,PATHWAY_type,by=\"#Term\")[,c(\"#Term\",\"ratio\",\"pvalue\",\"Corrected P-Value\",\"typeII\",\"typeI\")]
	if(nrow(PATHWAY_enrichraw)==0){
	        print(\"There is no enriched KEGG Pathways,please check your input file!\")
	}else{
		### color card:
		ramp_cc<-colorRamp(as.vector(c(\"#008B45\",\"white\")))
		continue_cols<-paste(rgb(ramp_cc(round(seq(0,1,length=10010),4)),max=255),\"E5\",sep=\"\")
		PATHWAY_enrich_raw<-PATHWAY_enrich_type[order(PATHWAY_enrich_type[,6],PATHWAY_enrich_type[,5],as.numeric(PATHWAY_enrich_type[,4]),as.numeric(PATHWAY_enrich_type[,3]),as.numeric(PATHWAY_enrich_type[,2])),]
                ### top bars to be view
		if(nrow(PATHWAY_enrich_raw)>60){
			topNs<-which(PATHWAY_enrich_raw[[3]]<0.5)            ## pvalue >0.5 were ignored 
				PATHWAY_enrich<-PATHWAY_enrich_raw[topNs,]
		}else{PATHWAY_enrich<-PATHWAY_enrich_raw}

		RichFactor<-round(PATHWAY_enrich[,2],2)            
		### significant stars:
		Pvalue<-as.numeric(PATHWAY_enrich[[3]]) ## pvalue
		#Pvalue<-as.numeric(PATHWAY_enrich[[4]])  ## fdr
		n_3tars<-which(Pvalue<=0.001)
		n_2tars<-setdiff(which(Pvalue<=0.01),n_3tars)
		n_1tars<-setdiff(which(Pvalue<=0.05),c(n_3tars,n_2tars))
		barcols<-getMyColor(targets=Pvalue)
		### pathway class:
		name2typeI<-PATHWAY_enrich[,c(1,6)]
		shortName<-gsub(\"Organismal Systems\",\"OS\",gsub(\"Human Diseases\",\"HD\",gsub(\"Genetic Information Processing\",\"GIP\",gsub(\"Environmental Information Processing\",\"EIP\",gsub(\"Drug Development\",\"DD\",gsub(\"Cellular Processes\",\"CP\",gsub(\"Metabolism\",\"M\",name2typeI[[2]])))))))
		xlabs<-paste(PATHWAY_enrich[[1]],shortName,sep=\":   \")
		typeII<-c(\"KEGG Pathway Class:\",\"EIP: Environmental Information Processing\",\"GIP: Genetic Information Processing\",\"CP : Cellular Processes\",\"OS : Organismal Systems\",\"DD : Drug Development\",\"HD : Human Diseases\",\"M  : Metabolism\")
		### ploting
		pdfname <-\"$opts{i}.pathway.pdf\"
		pdf(file=pdfname,width=w,height=h)
		skip<-max(RichFactor)/50
		
		par(fig=c(0,0.92,0,1),mar = c($opts{mb},$opts{ml},4,0.1)+0.1,new=F)
		bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab=\"EnrichmentRatio: Sample_number/Background_number\",cex.lab=1,font.lab=1,col.lab=\"black\",border=F)
		text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=$opts{pc},font=1)
		#text(bar,rep(-0.01,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=$opts{pc},font=1)
		
		if(length(n_3tars)!=0){text(bar[n_3tars],RichFactor[n_3tars]+skip,\"***\",cex=$opts{sc},col=\"black\")}
		if(length(n_2tars)!=0){text(bar[n_2tars],RichFactor[n_2tars]+skip,\"**\",cex=$opts{sc},col=\"black\")}
		if(length(n_1tars)!=0){text(bar[n_1tars],RichFactor[n_1tars]+skip,\"*\",cex=$opts{sc},col=\"black\")}
		legend(\"topleft\",typeII,cex=$opts{lc},bty=\"n\",inset=$opts{is})
		
		### legend for bar colors:
		par(fig=c(0.91,1,0,0.9),mar = c($opts{mb},0.1,4,4)+0.1,new=T)
		Bars<-barplot(rep(0.5,10000),horiz=T,col=continue_cols[10000:1],border=continue_cols[10000:1],space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
		text(rep(0.6,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
		text(0,length(continue_cols)+500,\"Pvalue\",adj=0,xpd=T)
		dev.off()
		
	}
	
}
if(type==\"DISEASE\"){
	print(\"Disease barplot\")
	PATHWAY_type<-read.delim(\"/mnt/lustre-2/users/wangyan/rna/database/KEGG/pathwayId_pathwayName_typeII_typeI\",sep=\"\\t\",header=T,check.names=F)
        PATHWAY_enrichraw<-read.delim(\"$opts{i}.disease\",header=T,sep=\"\\t\",check.names=F)
	PATHWAY_enrich_type<-merge(PATHWAY_enrichraw,PATHWAY_type,by=\"#Term\")[,c(\"#Term\",\"ratio\",\"pvalue\",\"Corrected P-Value\",\"typeII\",\"typeI\")]
	if(nrow(PATHWAY_enrichraw)==0){
	        print(\"There is no enriched KEGG Pathways,please check your input file!\")
	}else{
		### color card:
		ramp_cc<-colorRamp(as.vector(c(\"#008B45\",\"white\")))
		continue_cols<-paste(rgb(ramp_cc(round(seq(0,1,length=10010),4)),max=255),\"E5\",sep=\"\")
		PATHWAY_enrich_raw<-PATHWAY_enrich_type[order(PATHWAY_enrich_type[,6],PATHWAY_enrich_type[,5],as.numeric(PATHWAY_enrich_type[,4]),as.numeric(PATHWAY_enrich_type[,3]),as.numeric(PATHWAY_enrich_type[,2])),]
                ### top bars to be view
		if(nrow(PATHWAY_enrich_raw)>60){
			topNs<-which(PATHWAY_enrich_raw[[3]]<0.5) ## pvalue >0.5 were ignored 
			PATHWAY_enrich<-PATHWAY_enrich_raw[topNs,]
		}else{
			PATHWAY_enrich<-PATHWAY_enrich_raw
		}

		RichFactor<-round(PATHWAY_enrich[,2],2)            
		### significant stars:
		Pvalue<-as.numeric(PATHWAY_enrich[[3]]) ## pvalue
		#Pvalue<-as.numeric(PATHWAY_enrich[[4]])  ## fdr
		n_3tars<-which(Pvalue<=0.001)
		n_2tars<-setdiff(which(Pvalue<=0.01),n_3tars)
		n_1tars<-setdiff(which(Pvalue<=0.05),c(n_3tars,n_2tars))
		barcols<-getMyColor(targets=Pvalue)
		### pathway class:
		name2typeI<-PATHWAY_enrich[,c(1,6)]
		shortName<-gsub(\"Organismal Systems\",\"OS\",gsub(\"Human Diseases\",\"HD\",gsub(\"Genetic Information Processing\",\"GIP\",gsub(\"Environmental Information Processing\",\"EIP\",gsub(\"Drug Development\",\"DD\",gsub(\"Cellular Processes\",\"CP\",gsub(\"Metabolism\",\"M\",name2typeI[[2]])))))))
		xlabs<-paste(PATHWAY_enrich[[1]],shortName,sep=\":   \")
		typeII<-c(\"Pathway Class:\",\"EIP: Environmental Information Processing\",\"GIP: Genetic Information Processing\",\"CP : Cellular Processes\",\"OS : Organismal Systems\",\"DD : Drug Development\",\"HD : Human Diseases\",\"M  : Metabolism\")
		### ploting
		pdfname <-\"$opts{i}.disease.pdf\"
		pdf(file=pdfname,width=w,height=h)
		skip<-max(RichFactor)/50
		
		par(fig=c(0,0.92,0,1),mar = c($opts{mb},$opts{ml},4,0.1)+0.1,new=F)
		bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab=\"EnrichmentRatio: Sample_number/Background_number\",cex.lab=1,font.lab=1,col.lab=\"black\",border=F)
		#text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=$opts{pc},font=1)
		#text(bar,rep(-0.01,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=$opts{pc},font=1)
		
		if(length(n_3tars)!=0){text(bar[n_3tars],RichFactor[n_3tars]+skip,\"***\",cex=$opts{sc},col=\"black\")}
		if(length(n_2tars)!=0){text(bar[n_2tars],RichFactor[n_2tars]+skip,\"**\",cex=$opts{sc},col=\"black\")}
		if(length(n_1tars)!=0){text(bar[n_1tars],RichFactor[n_1tars]+skip,\"*\",cex=$opts{sc},col=\"black\")}
		legend(\"topleft\",typeII,cex=$opts{lc},bty=\"n\",inset=$opts{is})
		
		### legend for bar colors:
		par(fig=c(0.91,1,0,0.9),mar = c($opts{mb},0.1,4,4)+0.1,new=T)
		Bars<-barplot(rep(0.5,10000),horiz=T,col=continue_cols[10000:1],border=continue_cols[10000:1],space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
		text(rep(0.6,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
		text(0,length(continue_cols)+500,\"Pvalue\",adj=0,xpd=T)
		dev.off()
		
	}
}
if(type==\"GO\"){
        GO_enrichraw<-read.delim(\"$opts{i}.go.draw\",header=T,sep=\"\\t\",check.names=F)
	if(nrow(GO_enrichraw)==0){
	        print(\"There is no enriched GO Term,please check your input file!\")
	}else{
                GO_enrich<-GO_enrichraw[order(GO_enrichraw[,4],GO_enrichraw[,3],GO_enrichraw[,5],-GO_enrichraw[,2],decreasing=F),]
		RichFactor<-round(GO_enrich[,2],2)
		### color card:
		ramp_cc<-colorRamp(as.vector(c(\"#008B45\",\"white\")))
		continue_cols<-paste(rgb(ramp_cc(round(seq(0,1,length=10010),4)),max=255),\"E5\",sep=\"\")
		### significant stars:
		Pvalue<-as.numeric(GO_enrich[[5]])  ## pvalue
		#Pvalue<-as.numeric(GO_enrich[[3]])   ## fdr
		n_3tars<-which(Pvalue<=0.001)
		n_2tars<-setdiff(which(Pvalue<=0.01),n_3tars)
		n_1tars<-setdiff(which(Pvalue<=0.05),c(n_3tars,n_2tars))
		barcols<-getMyColor(targets=Pvalue)
		### pathway class:
		name2typeI<-GO_enrich[,c(1,4)]
		shortName<-gsub(\"biological_process\",\"BP\",gsub(\"cellular_component\",\"CC\",gsub(\"molecular_function\",\"MF\",name2typeI[[2]])))
		
		xlabs<-paste(GO_enrich[[1]],shortName,sep=\":   \")
		xlabsMax<-max(nchar(xlabs))
		scale<-xlabsMax/60.0
		typeII<-c(\"GO Class:\",\"BP: Biological Process\",\"CC: Cellular Component\",\"MF : Molecular Function\")
		### ploting
		pdfname <-\"$opts{i}.go.pdf\"
		pdf(file=pdfname,width=ceiling(w*scale),height=ceiling(h*scale))
		skip<-max(RichFactor)/50
		
		par(fig=c(0,0.92,0,1),mar = c(ceiling(mb*scale/1.1),ceiling(ml*scale/1.2),4,0.1)+0.1,new=F)
		bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab=\"EnrichmentRatio: Sample_number/Background_number\",cex.lab=1,font.lab=1,col.lab=\"black\",border=F)
		text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=$opts{pc},font=1)
		#text(bar,rep(-0.01,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=$opts{pc},font=1)
		
		if(length(n_3tars)!=0){text(bar[n_3tars],RichFactor[n_3tars]+skip,\"***\",cex=$opts{sc},col=\"black\")}
		if(length(n_2tars)!=0){text(bar[n_2tars],RichFactor[n_2tars]+skip,\"**\",cex=$opts{sc},col=\"black\")}
		if(length(n_1tars)!=0){text(bar[n_1tars],RichFactor[n_1tars]+skip,\"*\",cex=$opts{sc},col=\"black\")}
		
		
		legend(\"topleft\",typeII,cex=$opts{lc},bty=\"n\",inset=$opts{is})
		
		### legend for bar colors:
		par(fig=c(0.91,1,0,0.9),mar = c($opts{mb},0.1,4,4)+0.1,new=T)
		Bars<-barplot(rep(0.5,10000),horiz=T,col=continue_cols[10000:1],border=continue_cols[10000:1],space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
		#text(rep(0.6,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
		text(rep(0.6,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
		text(0,length(continue_cols)+500,\"Pvalue\",adj=0,xpd=T)
		#text(0,length(continue_cols)+500,\"FDR\",adj=0,xpd=T)
		dev.off()
		
	}
	## GO Plot End ##
}

";

#system ("R --restore --no-save < $opts{i}.r");
`R --restore --no-save < $opts{i}.r`;
#system ('rm *.r');
if(-e "$opts{i}.go"){system("rm $opts{i}.go")};
if(-e "$opts{i}.pathway"){system("rm $opts{i}.pathway")};
if(-e "$opts{i}.disease"){system("rm $opts{i}.disease")};
if(-e "$opts{i}.*.pathway.xls"){system("rm $opts{i}.*.pathway.xls")};
if(-e "$opts{i}.*.disease.xls"){system("rm $opts{i}.*.disease.xls")};

