#/usr/bin/env python
#coding:utf-8
#yanjun
#20170901

import argparse
import ConfigParser 
import os
import re
from PIL import Image,ImageDraw

parser=argparse.ArgumentParser(description="KEGG all protein/metabolic annotation")
parser.add_argument("-config",type=str,required=True,help="config.ini")
parser.add_argument("-pathwaytxt",type=str,required=True,help="contain acc2KO`s file")
parser.add_argument("-org",type=str,required=True,help='''
T01001|hsa|Homo sapiens|human|Eukaryotes|Animals|Vertebrates|Mammals;
T01710|gmx|Glycine max|soybean|Eukaryotes|Plants|Eudicots|Pea family;
T02846|ecol|Escherichia coli LY180|Prokaryotes|Bacteria|Gammaproteobacteria - Enterobacteria|Escherichia''')
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
		self.accs = self.get_acc()
		self.png_dir = png_dir
		self.html_dir = html_dir
		self.outdir = outdir
		
	def judge(self):
		'''判断这个通路是否有蛋白/代谢物注释上，返回注释蛋白个数'''
		mark = 0
		html_data = open(html_dir,'r').read()
		for acc in self.accs:
			if '%s ('%acc in html_data:
				'''每个acc都会有注释信息'''
				mark += 1
		return mark		
			
	def get_acc(self):
		'''将pathway.txt里面的对应信息找出来'''
		accs = [x.strip().split('\t')[1] for x in open(self.pathway,'r').readlines()]			
		return set(accs)
		
	def get_html(self):	
		'''绘制新的网页'''
		outdir = self.outdir
		name = self.name
		html_dir = self.html_dir 
		accs = self.accs
		with open('%s/%s.html'%(outdir,name),'w')as html_w:
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
					for acc in accs:
						if '%s ('%acc in line:
							new_line_1 = re.sub(r'href="/',r'href="http:/www.kegg.jp/',line)
							new_line_2 = re.sub(r' onmouseover=.*?/>',r'/>',new_line_1)
							html_w.write(new_line_2)
			html_w.write('''</map>
<div id="poplay" class="poplay" />
</form>
</body>
</html>''')			
			
		
	def get_png(self):
		'''对图片进行操作'''
		outdir = self.outdir
		name = self.name
		html_dir = self.html_dir 
		accs = self.accs	
		png = Image.open(self.png_dir)	
		draw = ImageDraw.Draw(png)
		html_data = open(self.html_dir,'r').read()
		for acc in accs:
			if '%s (' %acc in html_data:
				acc_data = re.findall(r'<area shape=(.*?)\tcoords=(.*?)\thref=.*?%s.*?" />'%acc,html_data)
				for i in range(len(acc_data)):
					shape = acc_data[i][0]
					coords = acc_data[i][1]
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
						draw.line([(x_left,y_left),(x_right,y_left),(x_right,y_right),(x_left,y_right),(x_left,y_left)],fill='red',width=1)
		png.save('%s/ko%s.png' %(outdir,name))
			
	def get_pathway_table(self):
		'''pathway_table.xls的单行输出'''
		html_data = open(self.html_dir,'r').read()
		DEFINITION = re.search(r'DEFINITION  (.*?)\n',html_data).group(1)
		acc_new = [acc for acc in self.accs if '%s (' %acc in html_data]
		url = 'http://www.kegg.jp/kegg-bin/show_pathway?ko%s' %(self.name)
		for acc in acc_new:
			url += '+%s' %(acc)
		line = 'ko%s\t%s\t%s\t%s\t%s\n' %(self.name,DEFINITION,len(acc_new),';'.join(acc_new),url)
		return line

if 1 == 2 :		
	with open(org_ko_list,'r')as org_r:
		for line in org_r.readlines():			
			pathway_ko = Pathway(line,pathwaytxt,)
			pathway_ko.get_png()
			pathway_ko.get_html()
			#print(pathway_ko.get_pathway_table())
	
if __name__ == '__main__':
	##python ~/develop_python_git/Annotation_GO.KEGG.COG/script/KEGG_annotation.py -config config.ini -pathwaytxt pathway.txt -org aaaaaa -out new
	pathway_00010 = Pathway('00010',pathwaytxt,'%s/ko00010.png'%pathwaysimage,'%s/ko00010.kgml'%pathwayshtml,args.out)
	print(pathway_00010.get_pathway_table())
	pathway_00010.get_html()
	print(u'网页生成成功')
	pathway_00010.get_png()
	print(u'图片生成成功')