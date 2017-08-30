#!/usr/bin/env python
import argparse

parser=argparse.ArgumentParser(description="get idmapping from large DB file used by mark")
parser.add_argument("-db",type=str,required=True,help="idmapping.dat|goa_uniprot_all.gaf")
parser.add_argument("-l1",type=int,required=True,help="acc1`s NR;for uniprot cols,-1 in py-script")
parser.add_argument("-l2",type=int,required=True,help="acc2`s NR;for GO|KO|COG cols ")
parser.add_argument("-mark",type=str,required=True,help="mark to specify line to used;GO:|KOG|KO")
parser.add_argument("-outfile",type=str,required=True,help="as uniprot2GO.list")

args=parser.parse_args()

uniprot2id = {}
out = open(args.outfile,'w')
with open(args.db,'r')as inf:
	for line in inf.readlines():
		if args.mark in line:
			uniprot = line.strip().split('\t')[args.l1-1]
			id = line.strip().split('\t')[args.l2-1]
			if not uniprot in uniprot2id.keys():
				for uniprot_old in uniprot2id.keys():
					out.write('%s\t%s\n' %(uniprot_old,uniprot2id[uniprot_old]))
					out.flush()
				uniprot2id = {}
				id_new = id
			else:
				id_old = uniprot2id[uniprot]
				id_new = '%s;%s' %(id_old,id)				
			uniprot2id[uniprot] = id_new

for uniprot_old in uniprot2id.keys():
	out.write('%s\t%s\n' %(uniprot_old,uniprot2id[uniprot_old]))
	out.flush()			

out.close()
