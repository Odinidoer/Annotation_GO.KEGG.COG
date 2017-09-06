#!/usr/bin/env python
#coding=utf-8
#20170906
#yanjun

import argparse
import ConfigParser 
import os
import re
import requests
import shutil
import time
import wget
from docx import Document
from multiprocessing import Pool
from random import choice


##获取当前时间和路径
now_time = time.strftime("%Y-%m-%d",time.localtime())
now_dir = os.getcwd()

##创建新的database文件夹
up_dir = os.path.split(now_dir)[0]
new_dir = '/'.join([up_dir,'database',now_time])

def os_mkdir(dir):
	if not os.path.isdir(dir):
		os.mkdir(dir)	
	
os_mkdir(new_dir)

##新建config.ini文件
with open('%s/config.ini' %(new_dir),'w') as config:
	config.write('''[Global]
last_update_time = {time}
author = yanjun

[idmapping]

goaother_url = ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz
idmapping_dat_url = ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
uniprot2GO.list = /mnt/ilustre/users/ting.kuang/scripts/proteomics/pipeline/bin/idmapping/uniprot2GO.list
uniprot2KO.list = /mnt/ilustre/users/ting.kuang/scripts/proteomics/pipeline/bin/idmapping/uniprot2ko.list
uniprot2COG.list = /mnt/ilustre/users/ting.kuang/scripts/proteomics/pipeline/bin/idmapping/uniprot2COG.list
uniprot2KOG.list = /mnt/ilustre/users/ting.kuang/scripts/proteomics/pipeline/bin/idmapping/uniprot2KOG.list

[GO]
go_obo = {dir}/go.obo

[KEGG]
#从 通路 到 蛋白/代谢物
ko1KO = {dir}/ko2KO.list

#从 蛋白/代谢物 到 通路
KO2ko = {dir}/KO2ko.list

#蛋白/代谢物的注释信息
KO2name = {dir}/KO2name.txt

#通路图片，ko背景
pathwaysimage = {dir}/pathwaysimage/

#通路具体信息，可用于网页显示
pathwayshtml = {dir}/pathwayshtml/

#KEGG富集，用于鉴定分类级别
class = {dir}/pathwayId_pathwayName_typeII_typeI

#pathway.txt 网页源代码，用于获取全部通路ID
pathway_all_html = {dir}/pathway.html

#用于获取KO的description 手动COPY全网页http://www.kegg.jp/kegg-bin/get_htext?ko00000.keg
KO_All_Categories = {dir}/KO-All-Categories.txt

#organism全物种的通路ko信息
organism2ko = {dir}/organism/

[COG]

[RELEASE]
release = {dir}/release.txt'''.format(time = now_time,dir = new_dir))

##获取新的配置文件信息
cfg = ConfigParser.ConfigParser()
cfg.read('%s/config.ini' %(new_dir))

##获取所有ko的列表
print('\npathway.list KO.list is downloading')
pathway_list_url = 'http://rest.kegg.jp/list/pathway'
pathway_list = re.findall(r'path:map(\d*)',requests.get(pathway_list_url).text)

##获取KEGG富集层次信息
pathway_ID_type_url = 'http://www.kegg.jp/kegg/pathway.html'
with open(cfg.get('KEGG','class'),'w')as class_out:
	pathway_data = requests.get(pathway_ID_type_url).text
	class_out.write('pathwayId\t#Term\ttypeII\ttypeI\n')
	for h4 in re.findall(r'<h4>\d.*? (.*?)</h4>(.*?)<hr class="frame3" />',pathway_data,re.S):
			for b in re.findall(r'<b>\d.*? (.*?)</b>(.*?)<div class="clear"></div>',h4[1],re.S):
					for a in re.findall(r'<dt>(\d*?)</dt><dd><.*?>(.*?)</a></dd>',b[1]):
							class_out.write('%s\t%s\t%s\t%s\n' %(a[0],a[1],b[0],h4[0]))
														
##获取所有KO/CID的列表的列表
COMPOUND_list_url = 'http://rest.kegg.jp/list/COMPOUND'
COMPOUND_data = requests.get(COMPOUND_list_url).text

ORTHOLOGY_list_url = 'http://rest.kegg.jp/list/ORTHOLOGY'
ORTHOLOGY_data = requests.get(ORTHOLOGY_list_url).text

with open(cfg.get("KEGG","KO2name"),'w')as KOname_w:
	KOname_w.write('%s\n%s' %(COMPOUND_data,ORTHOLOGY_data))

print('\npathway`s png and html is downloading')	
##下载所有ko的网页和png信息
pathwayshtml = cfg.get("KEGG","pathwayshtml")
os_mkdir(pathwayshtml)
pathwaysimage = cfg.get("KEGG","pathwaysimage")
os_mkdir(pathwaysimage)

if u'快速下载KEGG' == u'快速下载KEGG':
	download_ko_list = pathway_list
	def download_png_html(ko):
		ko = str(ko)
		try:
			kgml_url = 'http://www.kegg.jp/kegg-bin/show_pathway?ko%s' %ko	
			with open('%s/ko%s.kgml' %(pathwayshtml,ko),'w')as kgml:
				kgml.write(requests.get(kgml_url,timeout=120).text)	
			png_url = 'http://rest.kegg.jp/get/ko%s/image' %ko
			with open('%s/ko%s.png' %(pathwaysimage,ko),'wb')as png:
				png.write(requests.get(png_url,timeout=120).content)
		except:
			print('error',ko)
			return ko
			
	while len(download_ko_list)>0:
		print('download_ko_list remain %s' %len(download_ko_list))
		p = Pool(4)
		new_download_ko_list= p.map(download_png_html,download_ko_list)
		p.close()
		p.join()
		download_ko_list = [x for x in new_download_ko_list if x]
		print(download_ko_list)
	
print('\nall pathway`s png and html downloaded correctly')			
			
##获取所有物种信息
organism2ko = cfg.get("KEGG","organism2ko")
os_mkdir(organism2ko)
organism_list_url = 'http://rest.kegg.jp/list/organism'
all_organism_data = requests.get(organism_list_url).text
organism_list = [x.split('\t')[1] for x in all_organism_data.split('\n') if '\t' in x]
download_organism_list = organism_list	
	
if u'想快速下载物种' == u'想快速下载物种':
	def get_organism_ko_list(organism):
		organism = str(organism)
		try:			
			organism_ko_url = 'http://rest.kegg.jp/list/pathway/%s' %organism
			with open('%s/%s.txt' %(organism2ko,organism),'w')as org_file:
				org_file.write(requests.get(organism_ko_url,timeout=120).text)
		except:			
			print('error:',organism)
			return organism
	
	while len(download_organism_list)>0:
		print('download_organism_list remain %s' %len(download_organism_list))
		p = Pool(4)
		new_download_organism_list = p.map(get_organism_ko_list,download_organism_list)
		p.close()
		p.join()
		download_organism_list = [x for x in new_download_organism_list if x]
		print(download_organism_list)
		
print('\nall siginal org`s pathway.list downloaded correctly')			
		
##处理得到物种信息

def split_org_line(line):
	'''对每行进行分割，取每个元素'''
	org_list = []
	if '\t' in line:
		items = line.split('\t')
		org_list.append(str(items[1]))
		for x in str(items[3]).split(';'):
			org_list.append(x)
	return org_list

##得到全物种水平的所有元素	
org2family = []
for organism_item in all_organism_data.split('\n'):
	org2family += split_org_line(organism_item)
		
org2family = set(map(str,org2family))

def merge_org(org):
	'''对物种水平的物种元素进行数据汇总'''
	merge_list = []
	for organism_item in all_organism_data.split('\n'):
		org_list = split_org_line(organism_item)
		if org in org_list:
			merge_list.append(org_list[0])
	with open('%s/%s.merge.txt' %(organism2ko,org),'w') as org_file:
		for family in merge_list:
			data = open('%s/%s.txt' %(organism2ko,family),'r').read()
			org_file.write(data)
			
map(merge_org,org2family)	

def org_txt2list(org):
	'''将txt文件转换成list文件，中间进行去重复'''
	with open('%s/%s.list' %(organism2ko,org),'w') as org_file:
		txt_data = open('%s/%s.merge.txt' %(organism2ko,org),'r').read()
		data = re.findall(r'path:.*?(\d*?)\t',txt_data)
		org_file.write('\n'.join(set(map(str,data))))
				
map(org_txt2list,org2family)
cmd = 'cat %s/*.list |sort -u >%s/ALL.list' %(organism2ko,organism2ko)		
os.system(cmd)

###下载go.obo
print('\ngo.obo is downloading')
go_obo_url = 'http://purl.obolibrary.org/obo/go.obo'
go_obo_downloadname = wget.download(go_obo_url)
go_obo_name = cfg.get('GO','go_obo')
cmd = 'mv %s %s' %(go_obo_downloadname,go_obo_name)
os.system(cmd)
		
print('\n before docx')

##在结题报告中写入新的版本信息
with open(go_obo_name,'r')as go_obo_r:
	'''format-version: 1.2
	data-version: releases/2017-07-28'''
	format_version = go_obo_r.readline().strip()
	data_version = go_obo_r.readline().strip()
	
go_version = u'''此次项目中，GO数据更新日期是：%s,对应GO官方网站的数据格式是：%s,对应GO官方网站的数据更新日期是：%s。''' %(now_time,format_version,data_version)

kegg_release_url = 'http://www.kegg.jp/kegg/docs/relnote.html'
kegg_release_data = requests.get(kegg_release_url).text
kegg_release = re.findall(r'<h4>Current release</h4>\n\n(.*?)\n<ul>',kegg_release_data)
kegg_version = u'''此次项目中，KEGG数据更新日期是：%s,对应KEGG官方网站的数据格式是：%s'''%(now_time,kegg_release)

with open(cfg.get('RELEASE,release'),'w') as release_w:
	release_w.write(go_version+'\n'+kegg_version)

##链接config.ini
cmd = 'rm config.ini && ln -s %s/config.ini config.ini' %(new_dir)
os.system(cmd)