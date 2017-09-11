#!/usr/bin/env python
#coding:utf-8
#yanjun
#20170911

import argparse
import ConfigParser 
import sys
import re
from goatools.obo_parser import GODag

parser=argparse.ArgumentParser(description="GO all protein annotation")
parser.add_argument("-config",type=str,required=True,help="config.ini")
parser.add_argument("-i",type=str,required=True,help="GO_list")
parser.add_argument("-out",type=str,required=True,help="outdir GO_out.xls")
args=parser.parse_args()

cfg = ConfigParser.ConfigParser()
cfg.read(args.config)
obo_file = cfg.get("GO","go_obo")
go_table = args.i[2]
out_file = args.out[3]

go2parent = {}
go2name = {}
go2namespace = {}
go2def = {}
data = open(obo_file,'r').read()
Terms = data.split('[Term]')
All_lines = []
go2acc = {}
acc2go = {}
go2alt_id = {}

for term in Terms:
	if re.search(r'\nid: GO',term):
		GO_id = re.search(r'\nid: (GO:\d*)',term).groups()[0]
		name = re.search(r'name: (.*?)\n',term).groups()[0]
		go2name[GO_id] = name
		namespace = re.search(r'namespace: (.*?)\n',term).groups()[0]
		go2namespace[GO_id] = namespace 
		def_inf = re.search(r'def: (.*?)\n',term).groups()[0]
		go2def[GO_id] = def_inf
		if re.search(r'is_a: (.*?) ',term):
			go_parents = re.findall(r'is_a: (GO.*?) ',term)
			go2parent[GO_id] = ';'.join(go_parents)
		if re.search('alt_id',term):
			go_alts = re.findall('alt_id:(GO.*?)\n',term)
			for go_alt in go_alts:
				go2alt_id[go_alt] = GO_id
		if re.search(r'\nreplaced_by: ',term):
			go_alts = re.search('\nreplaced_by: (GO:\d*)',term).groups()[0]
			go2alt_id[GO_id] = go_alts

godata = open(go_table,'r').read()
GOs = re.findall(r'GO:\d*',godata)
GOS = set(GOs)

with open(go_table,'r')as gotable:
	for line in gotable.readlines():
		items = line.strip().split('\t')
		acc = items[0]
		gos = re.findall('GO:\d+',items[1])
		acc2go[acc] = gos
		for go in gos:
			if go in go2acc.keys():
				goacc = go2acc[go]
				goacc_new = goacc+';'+acc
			else:
				goacc_new = acc
			go2acc[go] = goacc_new	
		
for go in go2acc.keys():
	if go in go2alt_id.keys():
		go2acc[go2alt_id[go]] = go2acc[go]
				
def get_parents(lines):
	for line in lines:
		go = line.split('+')[-1]
		if go in go2parent.keys():
			parents = go2parent[go]
			for parent in parents.split(';'):
				line_new = '%s+%s' %(line,parent)
				lines.append(line_new)
			lines.remove(line)
	return lines		

for go in GOS:
	if go in go2alt_id.keys():
		go = go2alt_id[go]
	lines = ['%s'%go,]
	i = 1
	while i<20:
		i = i+1
		lines = get_parents(lines)
	All_lines += lines	

GO_all_dict_line = []

for line in All_lines:
	GO_all_dict_line += line.split('+')
GO_all_dict = set(GO_all_dict_line)

out = open(out_file,'w')
out.write('GO	level	depth	name	namespace	def	nums_of_GOs	GOs	num_of_Accs	Accs	Accs_GO\n')
GO_Parser = GODag(obo_file,optional_attrs=['relationship'])	
for i in range(1,20):
	for go in go2name.keys():
		go_raw = go
		if go in go2alt_id.keys() and go in GOS:
			go = go2alt_id[go_raw]
		if go in GO_all_dict:
			if GO_Parser[go].level is i-1:
				gos_set = set([line.split('+')[0] for line in All_lines if go in line])
				num_of_GO = len(gos_set)
				accs = ''
				for go_ac in gos_set:
					if accs == '':
						accs = go2acc[go_ac]
					else:	
						accs = accs+';'+go2acc[go_ac]
				accs = set(accs.split(';'))				
				num_of_Accs = len(accs)				
				accs_write = []				
				for acc in accs:				
					gos_from_acc = [gox for gox in acc2go[acc] if gox in gos_set]					
					accs_write.append('%s(%s)' %(acc,';'.join(gos_from_acc)))					
				out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' 
%(go_raw,GO_Parser[go].level+1,GO_Parser[go].depth+1,go2name[go],go2namespace[go],go2def[go],num_of_GO,';'.join(gos_set),num_of_Accs,';'.join(accs),';'.join(accs_write)))
				out.flush()	
				
out.close()