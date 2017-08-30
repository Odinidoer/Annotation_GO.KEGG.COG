#!/usr/bin/env python
#20170821
#yanjun
import argparse
import ConfigParser  
import requests
import re
import os

#config file
parser=argparse.ArgumentParser(description="config by config.ini")
parser.add_argument("-config",type=str,required=True,help="confirm config.ini")
args=parser.parse_args()  

cfg = ConfigParser.ConfigParser()
cfg.read(args.config)

pathwayshtml = cfg.get("KEGG","pathwayshtml")
KO_All_Categories = cfg.get("KEGG","KO_All_Categories")
ko2KO_out = open(cfg.get("KEGG","ko1KO"),'w')
KO2ko_out = open(cfg.get("KEGG","KO2ko"),'w')
class_out = open(cfg.get("KEGG","class"),'w')

ko_all = []
dict_KO2ko = {}
KO2desc = {}

with open(KO_All_Categories,'r')as KO_O:
	for line in KO_O.readlines():
		if ';' in line:
			items = line.split('  ')
			KO2desc[items[3].strip()] = items[4].strip()
			

pathway_data = open(cfg.get("KEGG","pathway_all_html"),'r').read()
class_out.write('pathwayId\t#Term\ttypeII\ttypeI\n')
for h4 in re.findall(r'<h4>\d.*? (.*?)</h4>(.*?)<hr class="frame3" />',pathway_data,re.S):
	for b in re.findall(r'<b>\d.*? (.*?)</b>(.*?)<div class="clear"></div>',h4[1],re.S):
		for a in re.findall(r'<dt>(\d*?)</dt><dd><.*?>(.*?)</a></dd>',b[1]):
			class_out.write('%s\t%s\t%s\t%s\n' %(a[0],a[1],b[0],h4[0]))
			ko_all.append(a[0])

def deal_ko_html(ko):
	KO_COs = []
	file = '%s/ko%s.kgml' %(pathwayshtml,ko)
	data = open(file).read()
	DEFINITION=re.search(r'DEFINITION  (.*)',data).group(1)
	
	#deal ONTOLOGE	
	if re.search(r'"/dbget-bin/www_bget',data):
		COM_INFO = re.findall(r'<area shape=circle.*?title="(.*?)"',data)
		for items in COM_INFO:
			for item in items.split(', '):
				CID = item[0:6]
				CIDdesc = item[8:-1]
				KO2desc[CID] = CIDdesc	
		KO_CO_info = re.findall(r'"/dbget-bin/www_bget.(.*?)"',data)
		for KO_CO_all in KO_CO_info:
			for KO_CO in KO_CO_all.split('+'):
				if 'C' in KO_CO or 'K' in KO_CO:	
					KO_COs.append(KO_CO)
					if KO_CO in dict_KO2ko.keys():
						val_new = dict_KO2ko[KO_CO]+';'+ko
					else:
						val_new = ko
					dict_KO2ko[KO_CO] = val_new
	ko2KO_out.write('%s\t%s\t%s\n' %(ko,DEFINITION,';'.join(set(KO_COs))))

map(deal_ko_html,ko_all)	

for KO_CO in dict_KO2ko.keys():
	if KO_CO in KO2desc.keys():
		KO2ko_out.write('%s\t%s\t%s\n' %(KO_CO,KO2desc[KO_CO],dict_KO2ko[KO_CO]))

ko2KO_out.close()	
KO2ko_out.close()
class_out.close()

##echo not mv
def sort_u_ko(X):
	cmd = 'mv %s %s.old && sort -u %s.old >%s && rm %s.old' %(X,X,X,X,X)
	os.system(cmd)	
sort_u_ko(cfg.get("KEGG","ko1KO"))	
sort_u_ko(cfg.get("KEGG","KO2ko"))	