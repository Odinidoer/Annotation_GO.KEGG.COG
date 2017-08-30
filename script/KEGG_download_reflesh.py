#!/usr/bin/env python
#20170815
#yanjun
import argparse
import ConfigParser  
import requests
import re

#config file
parser=argparse.ArgumentParser(description="config by config.ini")
parser.add_argument("-config",type=str,required=True,help="confirm config.ini")
parser.add_argument("-expt",type=str,help="expt.list")
args=parser.parse_args()  

cfg = ConfigParser.ConfigParser()
cfg.read(args.config)

pathwayshtml = cfg.get("KEGG","pathwayshtml")
pathwaysimage = cfg.get("KEGG","pathwaysimage")

pathways_all_url = 'http://www.kegg.jp/kegg/pathway.html'
kos_data = requests.get(pathways_all_url).text	
with open(cfg.get("KEGG","pathway_all_html"),'w')as pathway_file:
	pathway_file.write(kos_data)
kos = set(re.findall(r'<dt>(\d\d\d\d\d)</dt>',kos_data))

with open('%s/../all.kos.list' %pathwayshtml ,'w') as kos_list:
	kos_list.write('\n'.join(['ko%s' %ko for ko in kos]))
	
if args.expt:
	kos = open(args.expt,'r').read().split('\n')
	
def download_new_inf(ko):
	try:
		png_url = 'http://www.genome.jp/kegg/pathway/ko/ko%s.png' %ko
		with open('%s/ko%s.png' %(pathwaysimage,ko),'wb')as png:
			png.write(requests.get(png_url).content)
			
		kgml_url = 'http://www.kegg.jp/kegg-bin/show_pathway?map%s' %ko	
		with open('%s/ko%s.kgml' %(pathwayshtml,ko),'w')as kgml:
			kgml.write(requests.get(kgml_url).text)	
	except:
		print(ko)
	
map(download_new_inf,kos)	