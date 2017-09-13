#!/usr/bin/env python
#coding=utf-8
#20170906
#yanjun

import argparse
import ConfigParser
import sys
reload(sys)
sys.setdefaultencoding('utf-8')

parser = argparse.ArgumentParser(description="get cid from compound name in KEGG")
parser.add_argument("-config",type=str,required=True,help="config.ini")
parser.add_argument("-i",required=True,help="input file as meta.list")
parser.add_argument("-o",required=True,help="out file as name2cid.xls")
args = parser.parse_args()

dict_al2al = {
	u'α':u'alpha',
	u'β':u'beta',
	u'γ':u'gamma',
	u'δ':u'delta',
	u'ε':u'epsilon',
	u'ζ':u'zeta',
	u'η':u'eta',
	u'θ':u'theta',
	u'ι':u'iota',
	u'κ':u'kappa',
	u'λ':u'lambda',
	u'μ':u'mu',
	u'ν':u'nu',
	u'ξ':u'xi',
	u'ο':u'omicron',
	u'π':u'pi',
	u'ρ':u'rho',
	u'σ':u'sigma',
	u'τ':u'tau',
	u'υ':u'upsilon',
	u'φ':u'phi',
	u'χ':u'chi',
	u'ψ':u'psi',
	u'ω':u'omega'}

dic_name2cid = {}

cfg = ConfigParser.ConfigParser()
cfg.read(args.config)

KO2name = cfg.get('KEGG','KO2name')

def get_stardard(name):
	for alpha in dict_al2al.keys():
		if alpha in name:
			name = name.replace(alpha,dict_al2al[alpha])
	name_len = len(name)
	if name_len>0:
		name_head = name[0]
		name_tail = name[1:name_len+1]
		name = ''.join([name_head.upper(),name_tail])	
	return name
	
dic_name2cid = {}
with open(KO2name,'r')as KO2name_r:
	for line in KO2name_r.readlines():
		if '\t' in line:
			items = line.strip().split('\t')
			CID = items[0].split(':')[1]
			names = items[1].split('; ')
			for name in names:
				dic_name2cid[name] = CID
				dic_name2cid[get_stardard(name)] = CID
				
			
out_w = open(args.o,'w')

with open(args.i,'r')as in_r:
	for line in in_r.readlines():
		line = unicode(line.strip())
		if line in dic_name2cid.keys():
			out_w.write('%s\t%s\t%s\n' %(line,line,dic_name2cid[line]))
		elif get_stardard(line) in dic_name2cid.keys():
			out_w.write('%s\t%s\t%s\n' %(line,get_stardard(line),dic_name2cid[get_stardard(line)]))

out_w.close()		