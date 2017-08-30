#!/usr/bin/env python
#coding=utf-8
#20170830
#yanjun

import os
import time

#获取当前时间和路径
now_time = time.strftime("%Y-%m-%d-%H-%M-%S",time.localtime())
now_dir = os.getcwd()

#创建新的DB文件夹
up_dir = os.path.split(now_dir)[0]
new_dir = '/'.join([up_dir,'DB',now_time])
os.mkdir(new_dir)

def work_new():

	#新建config.ini文件
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

#通路图片，ko背景
pathwaysimage = {dir}/pathwaysimage/

#通路具体信息，可用于网页显示
pathwayshtml = {dir}/pathwayshtml/

#KEGG富集，用于鉴定分类级别
class = {dir}/pathwayId_pathwayName_typeII_typeI

#pathway.txt 网页源代码，用于获取全部通路ID
pathway_all_html = {dir}/pathway.html

#用于获取KO的description 手动COPY全网页http://www.kegg.jp/kegg-bin/get_htext?ko00000.keg
KO_All_Categories =  {dir}/KO-All-Categories.txt

[COG]
	'''.format(time = now_time,dir = new_dir))

	#链接config.ini
	cmd = 'rm %s/config.ini && ln -s %s/config.ini %s/config.ini' %(up_dir,new_dir,up_dir)
	os.system(cmd)
	
if __name__ == '__main__':
	work_new()