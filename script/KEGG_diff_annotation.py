#/usr/bin/env python
#coding:utf-8
#yanjun
#20170911

import argparse
import ConfigParser 
import os
import re
from PIL import Image,ImageDraw

parser=argparse.ArgumentParser(description="KEGG all protein/metabolic annotation")
parser.add_argument("-config",type=str,required=True,help="config.ini")
parser.add_argument("-pathwaytxt",type=str,required=True,help="contain acc2KO`s file")
parser.add_argument("-org",type=str,default='ALL',help='''hsa|Eukaryotes|Animals|Vertebrates|Mammals
gmx|Eukaryotes|Plants|Eudicots|Pea family
ecol|Prokaryotes|Bacteria|Gammaproteobacteria - Enterobacteria|Escherichia
default:ALL''')
parser.add_argument("-diffxls",type=str,required=True,help="as A_vs_B.diff.exp.xls")
parser.add_argument("-out",type=str,required=True,help="outdir for png|html|xls")

args=parser.parse_args()
cfg = ConfigParser.ConfigParser()
cfg.read(args.config)
organism2ko = cfg.get("KEGG","organism2ko")
org_ko_list = '%s/%s.list' %(organism2ko,args.org)
pathwayshtml = cfg.get("KEGG","pathwayshtml")
pathwaysimage = cfg.get("KEGG","pathwaysimage")

pathwaytxt = args.pathwaytxt
diffxls = args.diffxls

if not os.path.isdir(args.out):
	os.mkdir(args.out)

class Pathway_diff(object):
	'''差异蛋白注释分析'''
	
	def __init__(self,name,pathway,diff,png_dir,html_dir,outdir):
		self.name = name
		self.pathway = pathway
		self.diff = diff
		self.png_dir = png_dir
		self.html_dir = html_dir
		self.outdir = outdir
		self.KOs = self.get_KO()
		self.mark_KO = self.judge_KO()
		self.head = self.get_head()
		self.KO2acc = self.get_KO2acc()
		self.coord2KO = self.get_coord2KO()
		self.acc2inf = self.get_acc2inf()
		self.coord2color = self.get_coord2color()
		
	def get_KO(self):
		'''将pathway.txt里面的KO信息找出来'''
		KOs = [x.strip().split('\t')[1] for x in open(self.pathway,'r').readlines()]			
		return set(KOs)		
		
	def judge_KO(self):
		'''判断这个通路是否有蛋白/代谢物注释上，返回注释蛋白个数'''
		mark = 0
		html_data = open(self.html_dir,'r').read()
		for KO in self.KOs :
			if '%s ('%KO in html_data:
				'''每个KO都会有注释信息'''
				mark += 1
		return mark
		
	def get_head(self):
		'''判断这个通路是否有蛋白/代谢物注释上，返回注释蛋白个数'''
		head_data = open(self.diff,'r').readline()
		items = head_data.strip().split('\t')
		for i in range(len(items)):
			if 'regulate' in items[i]:
				mark_regulate = i
			if 'significant' in items[i]:
				mark_significant = i
			if 'FC' in items[i] and not 'log' in items[i]:
				mark_FC = i
		return(mark_regulate,mark_significant,mark_FC)
		
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
	
	def get_coord2KO(self):
		'''获取通路的全部可映射坐标，并取得坐标和蛋白/代谢的映射关系'''
		coord2KO = {}
		pathway_KOs = self.KOs
		html_data = open(self.html_dir,'r').read()
		coord2KO_info = re.findall(r'<area.*?coords=(.*?)\thref="/dbget-bin/www_bget.(.*?)"',html_data)
		for items in coord2KO_info:
			coord = items[0].replace(',','_')
			KOs_incoord = items[1]#join with '+'
			KOs_incoord = [KO for KO in KOs_incoord.split('+') if KO in pathway_KOs]
			if len(KOs_incoord) >0:
				coord2KO[coord] = '+'.join(KOs_incoord)
		return coord2KO
		
	def get_acc2inf(self,):
		'''通过差异表格得到上下调纤细颜色等信息'''
		acc2inf = {}
		with open(self.diff,'r')as diff_r:
			head = diff_r.readline()
			for line in diff_r.readlines():
				acc = line.split('\t')[0]
				inf = line.strip()
				acc2inf[acc] = inf
		return acc2inf
		
	def get_coord2color(self):
		'''获取每个可映射坐标的颜色，红色，绿色和混合色（半红半绿）'''
		KOs = self.KOs
		coord2KO = self.coord2KO
		KO2acc = self.KO2acc
		acc2inf = self.acc2inf
		mark = self.head
		colors = []
		coord2color = {}		
		for coord in coord2KO.keys():
			KO = coord2KO[coord].split('+')
			for one_KO in KO:
				acc = KO2acc[one_KO].split(';')
				for one_acc in acc:
					color_inf = acc2inf[one_acc]
					if 'yes' == color_inf.split('\t')[mark[1]]: 
						if 'up' == color_inf.split('\t')[mark[0]]:
							colors.append('red')
						if 'down' == color_inf.split('\t')[mark[0]]:
							colors.append('green')
			if len(colors)>0:				
				coord2color[coord] = ';'.join(colors)
		return coord2color
	
	def calc_html(self):
		'''得到结果文件网页'''		
		coord2KO = self.coord2KO
		coord2color = self.coord2color
		KO2acc = self.KO2acc
		mark = self.head
		FC_mark = int(mark[2])
		name = self.name
		acc2inf = self.acc2inf
		with open('%s/ko%s.html'%(self.outdir,name),'w')as html_w:
			html_w.write('''<html>
<head>
<meta http-equiv="content-type" content="text/html; charset=utf-8">
<title>ko%s.png</title>
<style type="text/css">
<!--
area {cursor: pointer;}
-->
</style>
<script type="text/javascript">
<!--
function showInfo(info) {
	obj = document.getElementById("result");
	obj.innerHTML = "<div style='cursor: pointer; position: absolute; right: 5px; color: #000;' onclick='javascript: document.getElementById(\"result\").style.display = \"none\";' title='close'>X</div>" + info;
	obj.style.top = document.body.scrollTop;
	obj.style.left = document.body.scrollLeft;
	obj.style.display = "";
}
//-->
</script>
</head>
<body>
<map name="ko%s.png">
'''%(name,name))
			for coord in coord2KO.keys():
				if coord in coord2color.keys():	
					if 'C' in coord2KO[coord]:
						html_w.write('''<area shape='circle' coords='%s' onmouseover='javascript: showInfo("<ul>'''%coord.replace('_',','))
						if re.subn('red','red',coord2color[coord])>0:
							html_w.write('''<li style=\"color: #f00;\">Up regulated	<ul>''')
							for KO in coord2KO[coord].split('+'):
								for acc in KO2acc[KO].split(';'):
									if 'up' in acc2inf[acc]:
										FC = acc2inf[acc].split('\t')[FC_mark]
										html_w.write('''<li>%s: %s (%s) </li>''' %(KO,acc,FC))
							html_w.write('''</ul></li>''')
						if re.subn('green','green',coord2color[coord])>0:
							html_w.write('''<li style=\"color: #0f0;\">Down regulated <ul>''')
							for KO in coord2KO[coord].split('+'):
								for acc in KO2acc[KO].split(';'):
									if 'down' in acc2inf[acc]:
										FC = acc2inf[acc].split('\t')[FC_mark]
										html_w.write('''<li>%s: %s (%s) </li>''' %(KO,acc,FC))
							html_w.write('''</ul></li>''')			
						html_w.write('''</ul>");' />''')			
					if 'K' in coord2KO[coord]:
						html_w.write('''<area shape='rect' coords='%s' onmouseover='javascript: showInfo("'''%(coord.replace('_',',')))
						if re.subn('red','red',coord2color[coord])>0:
							html_w.write('''<li style=\"color: #f00;\">Up regulated	<ul>''')
							for KO in coord2KO[coord].split('+'):
								for acc in KO2acc[KO].split(';'):
									if 'up' in acc2inf[acc]:
										FC = acc2inf[acc].split('\t')[FC_mark]
										html_w.write('''<li>%s: %s (%s) </li>''' %(KO,acc,FC))
							html_w.write('''</ul></li>''')		
						if re.subn('green','green',coord2color[coord])>0:
							html_w.write('''<li style=\"color: #0f0;\">Down regulated <ul>''')
							for KO in coord2KO[coord].split('+'):
								for acc in KO2acc[KO].split(';'):
									if 'down' in acc2inf[acc]:
										FC = acc2inf[acc].split('\t')[FC_mark]
										print(FC,acc2inf[acc])
										html_w.write('''<li>%s: %s (%s) </li>''' %(KO,acc,FC))
							html_w.write('''</ul></li>''')						
						html_w.write('''</ul>");' />''')			
			html_w.write('''
</map>
<img src='./ko%s.png' usemap='ko%s.png' />
<div id='result' style='position: absolute; width: 50%%; border: 1px solid #000; background-color: #fff; filter: alpha(opacity=95); opacity: 0.95; font-size: 12px; padding-right: 20px; display: none;' onmouseover="javascript: this.style.filter = 'alpha(opacity=100)'; this.style.opacity = 1;" onmouseout="javascript: this.style.filter ='alpha(opacity=95)'; this.style.opacity = 0.95;"></div>
</body></html>'''%(name,name))				
		
	def calc_png(self):
		'''得到结果文件图片'''
		coord2KO = self.coord2KO
		coord2color = self.coord2color
		png = Image.open(self.png_dir)	
		draw = ImageDraw.Draw(png)
		for coord in coord2KO.keys():
			items = coord.split('_')
			if 'C' in coord2KO[coord] and coord in coord2color.keys():	
				c_x = float(items[0])/1.
				c_y = float(items[1])/1.
				diameter = float(items[2])/1.
				x_left_up = c_x - diameter
				y_left_up = c_y - diameter
				x_right_down = c_x + diameter
				y_right_down = c_y + diameter
				if 'red' in coord2color[coord] and not 'green' in coord2color[coord]:	
					draw.ellipse((x_left_up,y_left_up,x_right_down,y_right_down),fill='red',outline='red')
				if 'green' in coord2color[coord] and not 'red' in coord2color[coord]:
					draw.ellipse((x_left_up,y_left_up,x_right_down,y_right_down),fill='green',outline='green')
				else:
					#draw.arc((x_left_up,y_left_up,x_right_down,y_right_down),0,180,fill='red') 
					#draw.arc((x_left_up,y_left_up,x_right_down,y_right_down),0,180,fill='green') 
					print('redgreenpoint')
			if 'K' in coord2KO[coord] and coord in coord2color.keys():
				x_left = float(items[0])
				y_left = float(items[1])
				x_right = float(items[2])
				y_right = float(items[3])
				if 'red' in coord2color[coord] and not 'green' in coord2color[coord]:
					draw.line([(x_left,y_left),(x_right,y_left),(x_right,y_right),(x_left,y_right),(x_left,y_left)],fill='red',width=2)
				if 'green' in coord2color[coord] and not 'red' in coord2color[coord]:
					draw.line([(x_left,y_left),(x_right,y_left),(x_right,y_right),(x_left,y_right),(x_left,y_left)],fill='green',width=2)
				else:
					draw.line([(x_left,y_left),(x_right,y_left),(x_right,y_right)],fill='red',width=2)
					draw.line([(x_right,y_right),(x_left,y_right),(x_left,y_left)],fill='green',width=2)
		png.save('%s/ko%s.png' %(self.outdir,self.name))
		
	def print_pathway_table(self):
		'''输出一行，pathway.table.xls'''
		pass
		
	def print_kegg_table(self):
		'''输出一行，kegg.table.xls'''
		pass

with open(org_ko_list,'r')as org_r:	
	for line in org_r.readlines():
		line = line.strip()
		Pathway = Pathway_diff(line,pathwaytxt,diffxls,'%s/ko%s.png'%(pathwaysimage,line),'%s/ko%s.kgml'%(pathwayshtml,line),args.out)
		if Pathway.mark_KO > 0:
			Pathway.calc_png()
			Pathway.calc_html()