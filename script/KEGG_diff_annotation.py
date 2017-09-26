#/usr/bin/env python
#coding:utf-8
#yanjun
#20170911

import argparse
import ConfigParser 
import os
import re
try:
	from PIL import Image,ImageDraw
except:
	import sys
	sys.path.append('/mnt/ilustre/users/jun.yan/.local/lib/python2.7/site-packages/Pillow-2.2.1-py2.7-linux-x86_64.egg')
	from PIL import Image,ImageDraw

parser=argparse.ArgumentParser(description="KEGG all protein/metabolic annotation")
parser.add_argument("-config",type=str,required=True,help="config.ini")
parser.add_argument("-pathwaytxt",type=str,required=True,help="contain acc2KO`s file")
parser.add_argument("-org",type=str,default='ALL-org',help='''hsa|Eukaryotes|Animals|Vertebrates|Mammals
gmx|Eukaryotes|Plants|Eudicots|Pea family
ecol|Prokaryotes|Bacteria|Gammaproteobacteria - Enterobacteria|Escherichia
default:ALL-org''')
parser.add_argument("-diffxls",type=str,required=True,help="as A_vs_B.diff.exp.xls")
parser.add_argument("-uplist",type=str,required=True,help="as A_vs_B.up.list")
parser.add_argument("-downlist",type=str,required=True,help="as A_vs_B.down.list")
parser.add_argument("-out",type=str,required=True,help="outdir for png|html|xls")

args=parser.parse_args()
cfg = ConfigParser.ConfigParser()
cfg.read(args.config)
organism2ko = cfg.get("KEGG","organism2ko")
org_ko_list = '%s/%s.list' %(organism2ko,args.org)
pathwayshtml = cfg.get("KEGG","pathwayshtml")
pathwaysimage = cfg.get("KEGG","pathwaysimage")

KO2acc = {}
acc2color = {}		
acc2inf = {}
#对三个字典开始赋值
with open(args.pathwaytxt,'r')as pathwaytxt_r:
	for line in pathwaytxt_r.readlines():
		items = line.strip().split('\t')
		acc = items[0]
		KO = items[1]
		if not KO in KO2acc.keys():
			new_accs = acc
		else:
			new_accs = '%s;%s' %(KO2acc[KO],acc)
		KO2acc[KO] = new_accs 

with open(args.uplist,'r')as uplist_r:
	for line in uplist_r.readlines():
		acc = line.strip()
		acc2color[acc] = 'red'
		
with open(args.downlist,'r')as downlist_r:
	for line in downlist_r.readlines():
		acc = line.strip()
		acc2color[acc] = 'green'
		
with open(args.diffxls,'r')as diffxls_r:
	head = diffxls_r.readline()
	head_items = head.strip().split('\t')
	for i in range(len(head_items)):
		if u'FC' in head_items[i] and not u'log' in head_items[i]:
			mark = i
	for line in diffxls_r.readlines():
		items = line.strip().split('\t')
		if mark < len(items)+1:
			acc = items[0]
			FC = items[mark]
			acc2inf[acc] = FC		
		
if not os.path.isdir(args.out):
	os.mkdir(args.out)
	
class Pathway_diff(object):
	'''差异蛋白KEGGt通路注释分析'''
	
	def __init__(self,name,KO2acc,acc2color,acc2inf,png_path,html_path,outdir):
		self.name = name
		self.KO2acc = KO2acc
		self.acc2color = acc2color
		self.acc2inf = acc2inf
		self.png_path = png_path
		self.html_path = html_path
		self.outdir = outdir
		self.coord2inf = self.get_coord2inf()
		self.coord2KO = self.coord2inf[0]
		self.coord2color = self.coord2inf[1]
		self.coord2acc = self.coord2inf[2]
		
	def get_coord2inf(self):
		'''得到某坐标的颜色信息，并用于返回参数判断此pathway是否被差异注释'''
		KO2acc = self.KO2acc
		acc2color = self.acc2color
		coord2KO = {}
		coord2acc = {}
		coord2color = {}
		html_data = open(self.html_path,'r').read()
		coord2KO_info = re.findall(r'<area.*?coords=(.*?)\thref="/dbget-bin/www_bget.(.*?)"',html_data)
		for items in coord2KO_info:
			KOs = []
			colors = []
			accs = []
			coord = items[0].replace(',','_')
			KOs_incoord = items[1]#join with '+'
			for KO in KOs_incoord.split('+'):
				if KO in KO2acc.keys():					
					KOs.append(KO)
					for acc in KO2acc[KO].split(';'):
						if acc in acc2color.keys():
							accs.append(acc)
							colors.append(acc2color[acc])
			if KOs:				
				coord2KO[coord] = ';'.join(KOs)
			if accs:
				coord2acc[coord] = ';'.join(accs)
			if colors:
				coord2color[coord] = ';'.join(colors)
		return(coord2KO,coord2color,coord2acc)
		
	def calc_html(self):
		'''得到结果文件网页'''	
		name = self.name
		KO2acc = self.KO2acc
		acc2color = self.acc2color
		acc2inf = self.acc2inf
		coord2KO = self.coord2KO
		coord2color = self.coord2color
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
	obj.innerHTML = "<div style='cursor: pointer; position: absolute; right: 5px; color: #000;' onclick='javascript: document.getElementById(\\"result\\").style.display = \\"none\\";' title='close'>X</div>" + info;
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
			for coord in coord2color.keys():
				up_num = re.subn('red','red',coord2color[coord])[1]
				down_num = re.subn('green','green',coord2color[coord])[1]
				if 'C' in coord2KO[coord]:
					html_w.write('''<area shape='circle' coords='%s' onmouseover='javascript: showInfo("<ul>'''%(coord.replace('_',',')))
					if up_num > 0:
						html_w.write('''<li style=\\"color: #f00;\\">Up regulated <ul>''')
						for KO in coord2KO[coord].split(';'):
							for acc in KO2acc[KO].split(';'):
								if acc in acc2color.keys():
									if 'red' == acc2color[acc]:
										FC = acc2inf[acc]
										html_w.write('''<li>%s: %s (%s) </li>''' %(KO,acc,FC))
						html_w.write('''</ul></li>''')		
					if down_num >0:
						html_w.write('''<li style=\\"color: #0f0;\\">Down regulated <ul>''')
						for KO in coord2KO[coord].split(';'):
							for acc in KO2acc[KO].split(';'):
								if acc in acc2color.keys():
									if 'green' in acc2color[acc]:
										FC = acc2inf[acc]
										html_w.write('''<li>%s: %s (%s) </li>''' %(KO,acc,FC))
						html_w.write('''</ul></li>''')						
					html_w.write('''</ul>");' />
''')
				if 'K' in coord2KO[coord]:
					html_w.write('''<area shape='rect' coords='%s' onmouseover='javascript: showInfo("<ul>'''%(coord.replace('_',',')))
					if up_num > 0:
						html_w.write('''<li style=\\"color: #f00;\\">Up regulated <ul>''')
						for KO in coord2KO[coord].split(';'):
							for acc in KO2acc[KO].split(';'):
								if acc in acc2color.keys():
									if 'red' == acc2color[acc]:
										FC = acc2inf[acc]
										html_w.write('''<li>%s: %s (%s) </li>''' %(KO,acc,FC))
						html_w.write('''</ul></li>''')		
					if down_num >0:
						html_w.write('''<li style=\\"color: #0f0;\\">Down regulated <ul>''')
						for KO in coord2KO[coord].split(';'):
							for acc in KO2acc[KO].split(';'):
								if acc in acc2color.keys():
									if 'green' in acc2color[acc]:
										FC = acc2inf[acc]
										html_w.write('''<li>%s: %s (%s) </li>''' %(KO,acc,FC))
						html_w.write('''</ul></li>''')						
					html_w.write('''</ul>");' />
''')			
			html_w.write('''
</map>
<img src='./ko%s.png' usemap='#ko%s.png' />
<div id='result' style='position: absolute; width: 50%%; border: 1px solid #000; background-color: #fff; filter: alpha(opacity=95); opacity: 0.95; font-size: 12px; padding-right: 20px; display: none;' onmouseover="javascript: this.style.filter = 'alpha(opacity=100)'; this.style.opacity = 1;" onmouseout="javascript: this.style.filter ='alpha(opacity=95)'; this.style.opacity = 0.95;"></div>
</body></html>'''%(name,name))				
		
	def calc_png(self):
		'''得到结果文件图片'''
		coord2KO = self.coord2KO
		coord2color = self.coord2color
		png = Image.open(self.png_path)	
		draw = ImageDraw.Draw(png)
		for coord in coord2KO.keys():
			items = coord.split('_')
			if coord in coord2color.keys():
				up_num = re.subn('red','red',coord2color[coord])[1]
				down_num = re.subn('green','green',coord2color[coord])[1]
				if 'C' in coord2KO[coord]: 
					c_x = float(items[0])/1.
					c_y = float(items[1])/1.
					diameter = float(items[2])/1.
					x_left_up = c_x - diameter
					y_left_up = c_y - diameter
					x_right_down = c_x + diameter
					y_right_down = c_y + diameter
					if up_num > 0 and down_num == 0:	
						draw.ellipse((x_left_up,y_left_up,x_right_down,y_right_down),fill='red',outline='red')
					elif up_num == 0 and down_num > 0:
						draw.ellipse((x_left_up,y_left_up,x_right_down,y_right_down),fill='green',outline='green')
					elif up_num > 0 and down_num > 0:
						draw.ellipse((x_left_up,y_left_up,x_right_down,y_right_down),fill='green',outline='blue')
						#draw.chord((x_left_up,y_left_up,x_right_down,y_right_down),0,180,fill='red') 
						#draw.chord((x_left_up,y_left_up,x_right_down,y_right_down),0,180,fill='green') 
				if 'K' in coord2KO[coord]:
					x_left = float(items[0])
					y_left = float(items[1])
					x_right = float(items[2])
					y_right = float(items[3])
					if up_num > 0 and down_num == 0:
						draw.line([(x_left,y_left),(x_right,y_left),(x_right,y_right),(x_left,y_right),(x_left,y_left)],fill='red',width=1)
					elif up_num == 0 and down_num > 0:
						draw.line([(x_left,y_left),(x_right,y_left),(x_right,y_right),(x_left,y_right),(x_left,y_left)],fill='green',width=1)
					elif up_num > 0 and down_num > 0:
						draw.line([(x_left,y_left),(x_right,y_left),(x_right,y_right)],fill='red',width=1)
						draw.line([(x_right,y_right),(x_left,y_right),(x_left,y_left)],fill='green',width=1)
		png.save('%s/ko%s.png' %(self.outdir,self.name))
		
	def print_table_line(self):
		new_kegg_table_line = []
		html_data = open(self.html_path,'r').read()
		DEFINITION = re.search(r'DEFINITION  (.*?)\n',html_data).group(1)
		coord2KO = self.coord2KO
		coord2acc = self.coord2acc
		acc2color = self.acc2color
		KO2acc = self.KO2acc
		url_raw = 'http://www.kegg.jp/kegg-bin/show_pathway?ko%s' %(self.name)
		accs_new = []
		KO_new = []
		for coord in coord2acc.keys():
			accs_new += coord2acc[coord].split(';')	
		accs_new = set(accs_new)
		for KO in KO2acc.keys():
			for acc in accs_new:
				if acc in KO2acc[KO]:
					KO_new.append(KO)
		KO_new = set(KO_new)
		KO2acc_str = []
		for KO in KO_new:
			KOsacc_diff = []
			for acc in KO2acc[KO].split(';'):
				if acc in accs_new:
					KOsacc_diff.append(acc)
			KO2acc_str.append('%s(%s)' %(KO,';'.join(KOsacc_diff)))
		KO2acc_str = ';'.join(KO2acc_str)
		pathway_url = url_raw
		for KO in KO_new:
			for acc in KO2acc[KO].split(';'):
				if acc in acc2color.keys():
					pathway_url += '+%s%%09%s' %(KO,acc2color[acc])			
		pathway_table_line = 'ko%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(self.name,DEFINITION,len(accs_new),';'.join(accs_new),len(KO_new),';'.join(KO_new),KO2acc_str,pathway_url)	
		for acc in accs_new:
			for KO in KO2acc.keys():
				if acc in KO2acc[KO]:
					kegg_url = '%s+%s%%09%s' %(url_raw,KO,acc2color[acc])
					if acc2color[acc] == 'red':
						up_down = 'up'
					elif acc2color[acc] == 'green':
						up_down = 'down'
					new_kegg_table_line.append('%s\t%s\t%s\tko%s\t%s\t%s\n' %(acc,up_down,KO,self.name,DEFINITION,kegg_url))
		kegg_table_line = ''.join(new_kegg_table_line)
		return(pathway_table_line,kegg_table_line)
		
	
pathway_table_w = open('%s/pathway_table.xls'%args.out,'w')
pathway_table_w.write('pathway\tpathway_name\tnumber_of_accs\taccs_list\tnumber_of_KOs\tKOs_list\tKO2acc\turl\n')
kegg_table_w = open('%s/kegg_table.xls'%args.out,'w')
kegg_table_w.write('acc\tregulate\tKO\tpathway\tpathway_name\turl\n')
with open(org_ko_list,'r')as org_r:	
	for line in org_r.readlines():
		line = line.strip()
		png_path = '%s/ko%s.png'%(pathwaysimage,line)
		html_path = '%s/ko%s.kgml'%(pathwayshtml,line)
		Pathway = Pathway_diff(line,KO2acc,acc2color,acc2inf,png_path,html_path,args.out)
		if Pathway.coord2color:
			Pathway.calc_png()
			Pathway.calc_html()	
			table_line = Pathway.print_table_line()
			pathway_table_w.write(table_line[0])
			kegg_table_w.write(table_line[1])
			
pathway_table_w.close()
kegg_table_w.close()			
