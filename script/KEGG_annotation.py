#/usr/bin/env python
#coding:utf-8
#yanjun
#20170906

import argparse
import ConfigParser 
import os
import re
from PIL import Image,ImageDraw

parser=argparse.ArgumentParser(description="KEGG all protein/metabolic annotation")
parser.add_argument("-config",type=str,required=True,help="config.ini")
parser.add_argument("-pathwaytxt",type=str,required=True,help="contain acc2KO`s file")
parser.add_argument("-org",type=str,required=True,help='''hsa|Eukaryotes|Animals|Vertebrates|Mammals
gmx|Eukaryotes|Plants|Eudicots|Pea family
ecol|Prokaryotes|Bacteria|Gammaproteobacteria - Enterobacteria|Escherichia''')
parser.add_argument("-out",type=str,required=True,help="outdir for png|html|xls")

args=parser.parse_args()
pathwaytxt = args.pathwaytxt
cfg = ConfigParser.ConfigParser()
cfg.read(args.config)
organism2ko = cfg.get("KEGG","organism2ko")
org_ko_list = '%s/%s.list' %(organism2ko,args.org)
pathwayshtml = cfg.get("KEGG","pathwayshtml")
pathwaysimage = cfg.get("KEGG","pathwaysimage")
if not os.path.isdir(args.out):
	os.mkdir(args.out)

class Pathway(object):
	'''对单个pathway进行操作，png/html,以及表格单行'''	
	
	def __init__(self,name,pathway,png_dir,html_dir,outdir):
		'''初始化类'''
		self.name = name
		self.pathway = pathway
		self.KOs = self.get_KO()
		self.KO2acc = self.get_KO2acc()
		self.png_dir = png_dir
		self.html_dir = html_dir
		self.outdir = outdir
		
	def judge(self):
		'''判断这个通路是否有蛋白/代谢物注释上，返回注释蛋白个数'''
		mark = 0
		html_data = open(self.html_dir,'r').read()
		for KO in self.KOs:
			if '%s ('%KO in html_data:
				'''每个KO都会有注释信息'''
				mark += 1
		return mark		
			
	def get_KO(self):
		'''将pathway.txt里面的KO信息找出来'''
		KOs = [x.strip().split('\t')[1] for x in open(self.pathway,'r').readlines()]			
		return set(KOs)
		
	def get_KO2acc(self):
		'''将pathway.txt里面的对应信息找出来'''
		KO2acc = {}
		for line in open(self.pathway,'r').readlines():
			items = line.strip().split('\t')
			if not items[1] in KO2acc.keys():
				KO2acc[items[1]] = items[0]
			else:
				KO2acc[items[1]] = KO2acc[items[1]]+';'+items[0]
		return KO2acc	
		
	def get_html(self):	
		'''绘制新的网页'''
		outdir = self.outdir
		name = self.name
		html_dir = self.html_dir 
		KOs = self.KOs
		with open('%s/ko%s.html'%(outdir,name),'w')as html_w:
			html_w.write('''<html>
<head>
<title>
ko%s.png
</title>
</head>
<body>
<img src="ko%s.png"  usemap="#mapdata" border="0" />
<map name="mapdata">
'''%(name,name))
			with open(html_dir,'r')as html_r:
				for line in html_r.readlines():
					for KO in KOs:
						if '%s ('%KO in line:
							new_line_1 = re.sub(r'href="/',r'href="http:/www.kegg.jp/',line)
							new_line_2 = re.sub(r' onmouseover=.*?/>',r'/>',new_line_1)
							html_w.write(new_line_2)
			html_w.write('''</map>
</body>
</html>''')			
			
		
	def get_png(self):
		'''对图片进行操作'''
		outdir = self.outdir
		name = self.name
		html_dir = self.html_dir 
		KOs = self.KOs	
		png = Image.open(self.png_dir)	
		draw = ImageDraw.Draw(png)
		html_data = open(self.html_dir,'r').read()
		for KO in KOs:
			if '%s (' %KO in html_data:
				KO_data = re.findall(r'<area shape=(.*?)\tcoords=(.*?)\thref=.*?%s.*?" />'%KO,html_data)
				for i in range(len(KO_data)):
					shape = KO_data[i][0]
					coords = KO_data[i][1]
					if shape == 'circle':
						items = coords.split(',')
						c_x = float(items[0])/1.
						c_y = float(items[1])/1.
						diameter = float(items[2])/1.
						x_left_up = c_x - diameter
						y_left_up = c_y - diameter
						x_right_down = c_x + diameter
						y_right_down = c_y + diameter
						draw.ellipse((x_left_up,y_left_up,x_right_down,y_right_down),fill='red',outline='red')
					if shape == 'rect':
						items = coords.split(',')
						x_left = float(items[0])
						y_left = float(items[1])
						x_right = float(items[2])
						y_right = float(items[3])
						draw.line([(x_left,y_left),(x_right,y_left),(x_right,y_right),(x_left,y_right),(x_left,y_left)],fill='red',width=2)
		png.save('%s/ko%s.png' %(outdir,name))
			
	def get_pathway_table(self):
		'''pathway_table.xls的单行输出'''
		html_data = open(self.html_dir,'r').read()
		DEFINITION = re.search(r'DEFINITION  (.*?)\n',html_data).group(1)
		KO_new = set([KO for KO in self.KOs if '%s (' %KO in html_data])
		accs_new = set(reduce((lambda x,y:x+'+'+y),[self.KO2acc[KO] for KO in KO_new]).split(';'))
		KO2acc_str = reduce((lambda x,y:x+';'+y),['%s(%s)' %(KO,self.KO2acc[KO]) for KO in KO_new])
		url = reduce((lambda x,y:x+'+'+y),KO_new,'http://www.kegg.jp/kegg-bin/show_pathway?ko%s' %(self.name))	
		line_pathway = 'ko%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(self.name,DEFINITION,len(accs_new),';'.join(accs_new),len(KO_new),';'.join(KO_new),KO2acc_str,url)
		return line_pathway

	def get_kegg_table(self):
		'''kegg_table.xls的多次单行输出'''
		html_data = open(self.html_dir,'r').read()
		DEFINITION = re.search(r'DEFINITION  (.*?)\n',html_data).group(1)
		KO_new = set([KO for KO in self.KOs if '%s (' %KO in html_data])
		line_kegg = ''
		for KO in KO_new:
			for acc in self.KO2acc[KO].split(';'):
				url = 'http://www.kegg.jp/kegg-bin/show_pathway?ko%s+%s' %(self.name,KO)
				line_kegg += '%s\t%s\tko%s\t%s\t%s\n' %(acc,KO,self.name,DEFINITION,url)
		return line_kegg		
		
			
		

pathway_table_w = open('%s/pathway_table.xls'%args.out,'w')
pathway_table_w.write('pathway\tpathway_name\tnumber_of_accs\taccs_list\tnumber_of_KOs\tKOs_list\tKO2acc\turl\n')
kegg_table_w = open('%s/kegg_table.xls'%args.out,'w')
kegg_table_w.write('acc\tKO\tpathway\tpathway_name\turl\n')		
with open(org_ko_list,'r')as org_r:	
	for line in org_r.readlines():
		line = line.strip()
		pathway_ko = Pathway(line,pathwaytxt,'%s/ko%s.png'%(pathwaysimage,line),'%s/ko%s.kgml'%(pathwayshtml,line),args.out)
		if pathway_ko.judge()>0:
			pathway_ko.get_png()
			pathway_ko.get_html()
			pathway_table_w.write(pathway_ko.get_pathway_table())
			kegg_table_w.write(pathway_ko.get_kegg_table())
			
pathway_table_w.close()
kegg_table_w.close()			