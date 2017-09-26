#/usr/bin/env python
#coding:utf-8
#yanjun
#20170913

import argparse
from scipy import stats

parser=argparse.ArgumentParser(description="enrichment diff pros")
parser.add_argument("-all",type=str,required=True,help="all pro`s KEGG pathway.table.xls")
parser.add_argument("-diff",type=str,required=True,help="diff pro`s KEGG pathway.table.xls")
parser.add_argument("-out_xls",type=str,required=True,help="outfile for enrichment.xls")
parser.add_argument("-out_plot",type=str,required=True,help="outfile for plot")

args=parser.parse_args()

ko2inf = {}
all_protein_list = []
diff_protein_list = []

with open(args.all,'r')as all_r:
	all_r.readline()
	'''pathway	pathway_name	number_of_accs	accs_list	number_of_KOs	KOs_list	KO2acc	url'''
	for line in all_r.readlines():
		items = line.strip().split('\t')
		ko = items[0]
		inf = line.strip()
		ko2inf[ko] = inf
		all_protein_list += items[3].split(';')
		
with open(args.diff,'r')as diff_r:
	diff_r.readline()
	for line in diff_r.readlines():
		items = line.strip().split('\t')
		diff_protein_list += items[3].split(';')
		
all_protein_list = set(all_protein_list)
diff_protein_list = set(diff_protein_list)

out_xls_w = open(args.out_xls,'w')
out_xls_w.write('pathway	pathway_name	all_number_of_accs	all_KO2acc	diff_number_of_accs	diff_KO2acc	p_value	url\n')
out_plot_w = open(args.out_plot,'w')
out_plot_w.write('''##Databases: KEGG PATHWAY								
##Statistical test method: hypergeometric test / Fisher's exact test								
##FDR correction method: Benjamini and Hochberg								
								
#Term	Database	ID	Input number	Background number	P-Value	Corrected P-Value	Input	Hyperlink\n''')

with open(args.diff,'r')as diff_r:
	diff_r.readline()
	'''pathway	pathway_name	number_of_accs	accs_list	number_of_KOs	KOs_list	KO2acc	url'''
	for line in diff_r.readlines():
		items = line.strip().split('\t')
		ko = items[0]
		ko_all_inf = ko2inf[ko].split('\t')
		all_number_of_accs = ko_all_inf[2]
		all_KO2acc = ko_all_inf[6]
		diff_number_of_accs = items[2]
		diff_KO2acc = items[6]
		x = int(diff_number_of_accs)#GO的差异蛋白
		m = int(all_number_of_accs)#GO的全部蛋白
		n = len(all_protein_list) - int(all_number_of_accs)#非GO的全部蛋白
		k = len(diff_protein_list)#次数
		p_value = stats.hypergeom.pmf(x,m+n,m,k)
		url = items[7]
		out_xls_w.write('%s\t%s\t%s/%s\t%s\t%s/%s\t%s\t%s\t%s\n' %(ko,items[1],all_number_of_accs,len(all_protein_list),all_KO2acc,diff_number_of_accs,len(diff_protein_list),diff_KO2acc,p_value,url))
		out_plot_w.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(items[1],'KEGG PATHWAY',ko,x,m,p_value,p_value,p_value,url))
		
out_plot_w.write('''--------------------								
								
#Term	Database	ID	Input number	Background number	P-Value	Corrected P-Value	Input	Hyperlink
								
--------------------								
								
#Term	Database	ID	Input number	Background number	P-Value	Corrected P-Value	Input	Hyperlink
								
--------------------''')		
		
out_xls_w.close()
out_plot_w.close()		
