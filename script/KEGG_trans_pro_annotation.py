#/usr/bin/env python
#coding:utf-8
#yanjun
#20171012

import argparse
import ConfigParser 
import os
import re
from collections import defaultdict
try:
	from PIL import Image,ImageDraw
except:
	import sys
	sys.path.append('/mnt/ilustre/users/jun.yan/.local/lib/python2.7/site-packages/Pillow-2.2.1-py2.7-linux-x86_64.egg')
	from PIL import Image,ImageDraw

parser=argparse.ArgumentParser(description="mRNA-pro KEGG pathway annotation")
parser.add_argument("-path",type=str,required=True,help="path dir saved png|html")
parser.add_argument("-mRNAs_pathway",type=str,required=True,help="mRNAs.pathway.txt")
parser.add_argument("-mRNAs_diff_xls",type=str,required=True,help="mRNAs.A_vs_B.diff.exp.xls")
parser.add_argument("-mRNAs_up",type=str,required=True,help="mRNAs.A_vs_B.up.list")
parser.add_argument("-mRNAs_down",type=str,required=True,help="mRNAs.A_vs_B.down.list")
parser.add_argument("-pros_pathway",type=str,required=True,help="mRNAs.pathway.txt")
parser.add_argument("-pros_diff_xls",type=str,required=True,help="pros.A_vs_B.diff.exp.xls")
parser.add_argument("-pros_up",type=str,required=True,help="pros.A_vs_B.up.list")
parser.add_argument("-pros_down",type=str,required=True,help="pros.A_vs_B.down.list")
parser.add_argument("-out",type=str,required=True,help="outdir for png|html|xls")
args = parser.parse_args()

#创建结果目录
if not os.path.isdir(args.out):
	pass
else:
	cmd = '/bin/rm -rf %s' %args.out
	os.system(cmd)
	
os.mkdir(args.out)
	
#创建字典存储数据
KO2rna = defaultdict(str)
rna2color = {}
rna2FC = {}
KO2pro = defaultdict(str)
pro2color = {}
pro2FC = {}

with open(args.mRNAs_pathway,'r')as pathwaytxt_r:
	for line in pathwaytxt_r.readlines():
		items = line.strip().split('\t')
		rna = items[0]
		KO = items[1]
		if not KO in KO2pro.keys():
			new_rnas = rna
		else:
			new_rnas = '%s;%s' %(KO2rna[KO],rna)
		KO2rna[KO] = new_rnas  

with open(args.mRNAs_up,'r')as uplist_r:
	for line in uplist_r.readlines():
		rna = line.strip()
		rna2color[rna] = 'red'
		
with open(args.mRNAs_down,'r')as downlist_r:
	for line in downlist_r.readlines():
		rna = line.strip()
		rna2color[rna] = 'green'
		
with open(args.mRNAs_diff_xls,'r')as diffxls_r:
	head = diffxls_r.readline()
	head_items = head.strip().split('\t')
	for i in range(len(head_items)):
		if u'FC' in head_items[i] and not u'log' in head_items[i]:
			mark = i
	for line in diffxls_r.readlines():
		items = line.strip().split('\t')
		if mark < len(items)+1:
			rna = items[0]
			FC = items[mark]
			rna2FC[rna] = FC
			if rna not in rna2color.keys():
				rna2color[rna] = 'yellow'				
		
with open(args.pros_pathway,'r')as pathwaytxt_r:
	for line in pathwaytxt_r.readlines():
		items = line.strip().split('\t')
		pro = items[0]
		KO = items[1]
		if not KO in KO2pro.keys():
			new_pros = pro
		else:
			new_pros = '%s;%s' %(KO2pro[KO],pro)
		KO2pro[KO] = new_pros		
		
with open(args.pros_up,'r')as uplist_r:
	for line in uplist_r.readlines():
		pro = line.strip()
		pro2color[pro] = 'red'
		
with open(args.pros_down,'r')as downlist_r:
	for line in downlist_r.readlines():
		pro = line.strip()
		pro2color[pro] = 'green'
		
with open(args.pros_diff_xls,'r')as diffxls_r:
	head = diffxls_r.readline()
	head_items = head.strip().split('\t')
	for i in range(len(head_items)):
		if u'FC' in head_items[i] and not u'log' in head_items[i]:
			mark = i
	for line in diffxls_r.readlines():
		items = line.strip().split('\t')
		if mark < len(items)+1:
			pro = items[0]
			FC = items[mark]
			pro2FC[pro] = FC
			if pro not in pro2color.keys():
				pro2color[pro] = 'yellow'	
				
def MRNA_pro_KEGG(name,indir,outdir,KO2rna,rna2color,rna2FC,KO2pro,pro2color,pro2FC):
	'''处理每一个通路'''
	path_combine = ''
	png_path = '%s/%s.png' %(indir,name)
	html_path = '%s/%s.html' %(indir,name)
	mark = 0
	try:
		background = Image.open(png_path)
	except:
		mark = 1
	if mark == 1:
		print ('%s.png or html is empty or not found!' %name)
		return 0
	
	KOs_set = ''
	coord_draw = Image.new("RGBA", background.size)
	layer = ImageDraw.Draw(coord_draw)
	html_data = open(html_path,'r').read()
	coord2KO_info = re.findall(r'<area shape=(.*?)	coords=(.*?)\thref="/dbget-bin/www_bget.(.*?)"',html_data)
	map_size = float(re.search(r'selected>(\d*?)\%',html_data).group(1))
	
	html_w = open('%s/%s.html' %(outdir,name),'w')
	html_w.write('''<html>
<head>
<meta http-equiv="content-type" content="text/html; charset=utf-8">
<title>%s.png</title>
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
<map name="%s.png">
'''%(name,name))
	

	
	for items_coord in coord2KO_info:
		rna_colors = []
		pro_colors = []
		
		shape = items_coord[0]
		items = items_coord[1].split(',')
		KOs_coord = items_coord[2].split('+')
		KOs = '+'.join(set([x for x in KO2pro.keys() if x in KOs_coord ] + [x for x in KO2rna.keys() if x in KOs_coord ]))
		KOs_set = '%s+%s' %(KOs_set,KOs)
		for KO in KOs.split('+'):
						
			'''左边转录颜色采集'''	
			if KO in KO2rna.keys():
				for rna in KO2rna[KO].split(';'):
					if rna in rna2color.keys():
						rna_colors.append(rna2color[rna])
					
			'''右边转录颜色采集'''
			if KO in KO2pro.keys():
				for pro in KO2pro[KO].split(';'):
					if pro in pro2color.keys():
						pro_colors.append(pro2color[pro])
				
		if len(rna_colors) or len(pro_colors):
			if shape == 'rect':
				'''矩形小方框'''
				html_w.write('''<area shape='rect' coords='%s' onmouseover='javascript: showInfo("<ul>'''%(','.join(items)))
				
				x_left = float(items[0])*100.0/map_size
				y_down = float(items[1])*100.0/map_size
				x_right = float(items[2])*100.0/map_size
				y_up = float(items[3])*100.0/map_size
				
				layer.rectangle((x_left+1.0, y_up+1.0, x_right-1.0, y_up-1.0), fill=(255, 255, 255, 128))
				
				if len(rna_colors):
					up_num = re.subn('red','red',';'.join(rna_colors))[1] #255,0,0
					down_num = re.subn('green','green',';'.join(rna_colors))[1] #0,255,0
					no_num = re.subn('yellow','yellow',';'.join(rna_colors))[1] #255,255,0
					
					if up_num > 0:
						html_w.write('''<li style=\\"color: #f00;\\">RNA up regulated <ul>''')
						for KO in KOs.split('+'):
							for rna in KO2rna[KO].split(';'):
								#print KO2rna[KO]
								if rna in rna2color.keys():
									if 'red' == rna2color[rna]:
										FC = rna2FC[rna]
										html_w.write('''<li>%s: %s (%s) </li>''' %(KO,rna,FC))
						html_w.write('''</ul></li>''')		
					if down_num >0:
						html_w.write('''<li style=\\"color: #006400;\\">RNA down regulated <ul>''')
						for KO in KOs.split('+'):
							for rna in KO2rna[KO].split(';'):
								if rna in rna2color.keys():
									if 'green' in rna2color[rna]:
										FC = rna2FC[rna]
										html_w.write('''<li>%s: %s (%s) </li>''' %(KO,rna,FC))
						html_w.write('''</ul></li>''')
					if no_num > 0:
						html_w.write('''<li style=\\"color: #CD6600;\\">RNA no regulated <ul>''')
						for KO in KOs.split('+'):
							for rna in KO2rna[KO].split(';'):
								if rna in rna2color.keys():
									if 'yellow' == rna2color[rna]:
										FC = rna2FC[rna]
										html_w.write('''<li>%s: %s (%s) </li>''' %(KO,rna,FC))
						html_w.write('''</ul></li>''')
					
					if up_num >0 and down_num == 0:
						'''只有上调和无明显表达'''
						left_x = x_left
						up_y = y_up
						right_x = x_left + (x_right -x_left)/2
						down_y = y_down
						layer.rectangle((left_x, up_y, right_x, down_y), fill=(255, 0, 0, 128))
						
					if up_num == 0 and down_num > 0:
						'''只有下调和无明显表达'''
						left_x = x_left
						up_y = y_up
						right_x = x_left +(x_right -x_left)/2
						down_y = y_down
						layer.rectangle((left_x, up_y, right_x, down_y), fill=(0, 255, 0, 128))	
						
					if up_num == 0 and down_num == 0 and no_num > 0:
						'''只有无明显表达'''
						left_x = x_left
						up_y = y_up
						right_x = x_left + (x_right -x_left)/2
						down_y = y_down
						layer.rectangle((left_x, up_y, right_x, down_y), fill=(255, 255, 0, 128))	
						
					if up_num > 0 and down_num > 0 and no_num == 0:
						'''只有上调和下调'''
						left_x = x_left
						up_y = y_up
						right_x = x_left + (x_right -x_left)/2
						down_y = (y_up + y_down )/2
						layer.rectangle((left_x, up_y, right_x, down_y), fill=(0, 255, 0, 128))	
						
						left_x = x_left
						up_y = (y_up + y_down )/2
						right_x = x_left + (x_right -x_left)/2
						down_y = y_down
						layer.rectangle((left_x, up_y, right_x, down_y), fill=(255, 0, 0, 128))
						
					if up_num > 0 and down_num > 0 and no_num > 0:
						'''上调下调无明显表达都有'''
						left_x = x_left
						up_y = y_up
						right_x = x_left + (x_right -x_left)/2
						down_y = y_up + (y_down - y_up)/3.0*1
						layer.rectangle((left_x, up_y, right_x, down_y), fill=(0, 255, 0, 128))	
						
						left_x = x_left
						up_y = y_up + (y_down - y_up)/3.0*1
						right_x = x_left + (x_right -x_left)/2
						down_y = y_up + (y_down - y_up)/3.0*2
						layer.rectangle((left_x, up_y, right_x, down_y), fill=(255, 255, 0, 128))
						
						left_x = x_left
						up_y = y_up + (y_down - y_up)/3.0*2
						right_x = x_left + (x_right -x_left)/2
						down_y = y_down
						layer.rectangle((left_x, up_y, right_x, down_y), fill=(255, 0, 0, 128))
						
				if len(pro_colors):
					up_num = re.subn('red','red',';'.join(pro_colors))[1] #255,0,0
					down_num = re.subn('green','green',';'.join(pro_colors))[1] #0,255,0
					no_num = re.subn('yellow','yellow',';'.join(pro_colors))[1] #255,255,0
					
					if up_num > 0:
						html_w.write('''<li style=\\"color: #f00;\\">Protein up regulated <ul>''')
						for KO in KOs.split('+'):
							for pro in KO2pro[KO].split(';'):
								if pro in rna2color.keys():
									if 'red' == pro2color[pro]:
										FC = pro2FC[pro]
										html_w.write('''<li>%s: %s (%s) </li>''' %(KO,pro,FC))
						html_w.write('''</ul></li>''')		
					if down_num >0:
						html_w.write('''<li style=\\"color: #006400;\\">Protein down regulated <ul>''')
						for KO in KOs.split('+'):
							for pro in KO2pro[KO].split(';'):
								if pro in pro2color.keys():
									if 'green' in pro2color[pro]:
										FC = pro2FC[pro]
										html_w.write('''<li>%s: %s (%s) </li>''' %(KO,pro,FC))
						html_w.write('''</ul></li>''')
					if no_num > 0:
						html_w.write('''<li style=\\"color: #CD6600;\\">Protein no regulated <ul>''')
						for KO in KOs.split('+'):
							for pro in KO2pro[KO].split(';'):
								if pro in pro2color.keys():
									if 'yellow' == pro2color[pro]:
										FC = pro2FC[pro]
										html_w.write('''<li>%s: %s (%s) </li>''' %(KO,pro,FC))
						html_w.write('''</ul></li>''')
					
					if up_num >0 and down_num == 0:
						'''只有上调和无明显表达'''
						left_x = x_left + (x_right -x_left)/2
						up_y = y_up
						right_x = x_right
						down_y = y_down
						layer.rectangle((left_x, up_y, right_x, down_y), fill=(255, 0, 0, 128))
						
					if up_num == 0 and down_num > 0:
						'''只有下调和无明显表达'''
						left_x = x_left + (x_right -x_left)/2
						up_y = y_up
						right_x = x_right
						down_y = y_down
						layer.rectangle((left_x, up_y, right_x, down_y), fill=(0, 255, 0, 128))	
						
					if up_num == 0 and down_num == 0 and no_num > 0:
						'''只有无明显表达'''
						left_x = x_left + (x_right -x_left)/2
						up_y = y_up
						right_x = x_right
						down_y = y_down
						layer.rectangle((left_x, up_y, right_x, down_y), fill=(255, 255, 0, 128))	
						
					if up_num > 0 and down_num > 0 and no_num == 0:
						'''只有上调和下调'''
						left_x = x_left + (x_right -x_left)/2
						up_y = y_up
						right_x = x_right
						down_y = (y_up + y_down )/2
						layer.rectangle((left_x, up_y, right_x, down_y), fill=(0, 255, 0, 128))	
						
						left_x = x_left + (x_right -x_left)/2
						up_y = y_up + (y_down -y_up)/2
						right_x = x_right
						down_y = y_down 
						layer.rectangle((left_x, up_y, right_x, down_y), fill=(255, 0, 0, 128))
						
					if up_num > 0 and down_num > 0 and no_num > 0:
						'''上调下调无明显表达都有'''
						left_x = x_left + (x_right -x_left)/2
						up_y = y_up
						right_x = x_right
						down_y = y_down + (y_up - y_down)/3.0*1
						layer.rectangle((left_x, up_y, right_x, down_y), fill=(0, 255, 0, 128))	
						
						left_x = x_left + (x_right -x_left)/2
						up_y = y_down + (y_up - y_down)/3.0*2
						right_x = x_right
						down_y = y_down + (y_up - y_down)/3.0*1
						layer.rectangle((left_x, up_y, right_x, down_y), fill=(255, 255, 0, 128))
						
						left_x = x_left + (x_right -x_left)/2
						up_y = y_down + (y_up - y_down)/3.0*1
						right_x = x_right
						down_y = y_down
						layer.rectangle((left_x, up_y, right_x, down_y), fill=(255, 0, 0, 128))
						
			if shape == 'poly':
				'''多边形，目前判断是长线条,有得是水平的，有的是长条的'''
				html_w.write('''<area shape='poly' coords='%s' onmouseover='javascript: showInfo("<ul>'''%(','.join(items)))
				
				x1 = float(items[0])*100.0/map_size
				y1 = float(items[1])*100.0/map_size
				x2 = float(items[2])*100.0/map_size
				y2 = float(items[3])*100.0/map_size
				x3 = float(items[4])*100.0/map_size
				y3 = float(items[5])*100.0/map_size
				x4 = float(items[6])*100.0/map_size
				y4 = float(items[7])*100.0/map_size
					
				layer.polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)], fill=(255, 255, 255, 128))	
				
				if len(rna_colors):
					up_num = re.subn('red','red',';'.join(rna_colors))[1] #255,0,0
					down_num = re.subn('green','green',';'.join(rna_colors))[1] #0,255,0
					no_num = re.subn('yellow','yellow',';'.join(rna_colors))[1] #255,255,0
					
					if up_num > 0:
						html_w.write('''<li style=\\"color: #f00;\\">RNA up regulated <ul>''')
						for KO in KOs.split('+'):
							for rna in KO2rna[KO].split(';'):
								if rna in rna2color.keys():
									if 'red' == rna2color[rna]:
										FC = rna2FC[rna]
										html_w.write('''<li>%s: %s (%s) </li>''' %(KO,rna,FC))
						html_w.write('''</ul></li>''')		
					if down_num >0:
						html_w.write('''<li style=\\"color: #006400;\\">RNA down regulated <ul>''')
						for KO in KOs.split('+'):
							for rna in KO2rna[KO].split(';'):
								if rna in rna2color.keys():
									if 'green' in rna2color[rna]:
										FC = rna2FC[rna]
										html_w.write('''<li>%s: %s (%s) </li>''' %(KO,rna,FC))
						html_w.write('''</ul></li>''')
					if no_num > 0:
						html_w.write('''<li style=\\"color: #CD6600;\\">RNA no regulated <ul>''')
						for KO in KOs.split('+'):
							for rna in KO2rna[KO].split(';'):
								if rna in rna2color.keys():
									if 'yellow' == rna2color[rna]:
										FC = rna2FC[rna]
										html_w.write('''<li>%s: %s (%s) </li>''' %(KO,rna,FC))
						html_w.write('''</ul></li>''')
					
					if up_num >0 and down_num == 0:
						'''只有上调和无明显表达'''
						x1_p = x1
						y1_p = y1
						x2_p = x2 + (x2 - x1)/2
						y2_p = y2 + (y2 - y1)/2
						x3_p = x3 + (x3 - x2)/2
						y3_p = y3 + (y3 - y2)/2
						x4_p = x4
						y4_p = y4
						
						layer.polygon([(x1_p,y1_p),(x2_p,y2_p),(x3_p,y3_p),(x4_p,y4_p)], fill=(255, 0, 0, 128))
						
					if up_num == 0 and down_num > 0:
						'''只有下调和无明显表达'''
						x1_p = x1
						y1_p = y1
						x2_p = x2 + (x2 - x1)/2
						y2_p = y2 + (y2 - y1)/2
						x3_p = x3 + (x3 - x2)/2
						y3_p = y3 + (y3 - y2)/2
						x4_p = x4
						y4_p = y4
						
						layer.polygon([(x1_p,y1_p),(x2_p,y2_p),(x3_p,y3_p),(x4_p,y4_p)], fill=(0, 255, 0, 128))
						
					if up_num == 0 and down_num == 0 and no_num > 0:
						'''只有无明显表达'''
						x1_p = x1
						y1_p = y1
						x2_p = x2 + (x2 - x1)/2
						y2_p = y2 + (y2 - y1)/2
						x3_p = x3 + (x3 - x2)/2
						y3_p = y3 + (y3 - y2)/2
						x4_p = x4
						y4_p = y4
						
						layer.polygon([(x1_p,y1_p),(x2_p,y2_p),(x3_p,y3_p),(x4_p,y4_p)], fill=(255, 255, 0, 128))
						
					if up_num > 0 and down_num > 0 and no_num == 0:
						'''只有上调和下调'''
						x1_p = x1
						y1_p = y1
						x2_p = x1 + (x2 - x1)/2
						y2_p = y1 + (y2 - y1)/2
						x5_p = x4
						y5_p = y4
						x6_p = x3 + (x4 - x3)/2
						y6_p = y3 + (y4 - y3)/2
						x3_p = x1_p + (x5_p - x1_p)/2
						y3_p = y1_p + (y5_p - y1_p)/2
						x4_p = x2_p + (x6_p - x2_p)/2
						y4_p = y2_p + (y6_p - y2_p)/2
						
						layer.polygon([(x1_p,y1_p),(x2_p,y2_p),(x4_p,y4_p),(x3_p,y3_p)], fill=(255, 0, 0, 128))
						layer.polygon([(x3_p,y3_p),(x4_p,y4_p),(x6_p,y6_p),(x5_p,y5_p)], fill=(0, 255, 0, 128))
						
					if up_num > 0 and down_num > 0 and no_num > 0:
						'''上调下调无明显表达都有'''
						x1_p = x1
						y1_p = y1
						x2_p = x1 + (x2 - x1)/2
						y2_p = y1 + (y2 - y1)/2
						x7_p = x4
						y7_p = y4
						x8_p = x3 + (x4 - x3)/2
						y8_p = y3 + (y4 - y3)/2
						x3_p = x1_p + (x7_p - x1_p)/3*1
						y3_p = y1_p + (y7_p - y1_p)/3*1
						x4_p = x2_p + (x8_p - x2_p)/3*1
						y4_p = y2_p + (y8_p - y2_p)/3*1
						x5_p = x1_p + (x7_p - x1_p)/3*2
						y5_p = y1_p + (y7_p - y1_p)/3*2
						x6_p = x2_p + (x8_p - x2_p)/3*2
						y6_p = y2_p + (y8_p - y2_p)/3*2
						
						layer.polygon([(x1_p,y1_p),(x2_p,y2_p),(x4_p,y4_p),(x3_p,y3_p)], fill=(255, 0, 0, 128))
						layer.polygon([(x3_p,y3_p),(x4_p,y4_p),(x6_p,y6_p),(x5_p,y5_p)], fill=(255, 255, 0, 128))
						layer.polygon([(x5_p,y5_p),(x6_p,y6_p),(x8_p,y8_p),(x7_p,y7_p)], fill=(0, 255, 0, 128))
						
				if len(pro_colors):
					up_num = re.subn('red','red',';'.join(pro_colors))[1] #255,0,0
					down_num = re.subn('green','green',';'.join(pro_colors))[1] #0,255,0
					no_num = re.subn('yellow','yellow',';'.join(pro_colors))[1] #255,255,0
					
					if up_num > 0:
						html_w.write('''<li style=\\"color: #f00;\\">Protein up regulated <ul>''')
						for KO in KOs.split('+'):
							for pro in KO2pro[KO].split(';'):
								if pro in rna2color.keys():
									if 'red' == pro2color[pro]:
										FC = pro2FC[pro]
										html_w.write('''<li>%s: %s (%s) </li>''' %(KO,pro,FC))
						html_w.write('''</ul></li>''')		
					if down_num >0:
						html_w.write('''<li style=\\"color: #006400;\\">Protein down regulated <ul>''')
						for KO in KOs.split('+'):
							for pro in KO2pro[KO].split(';'):
								if pro in pro2color.keys():
									if 'green' in pro2color[pro]:
										FC = pro2FC[pro]
										html_w.write('''<li>%s: %s (%s) </li>''' %(KO,pro,FC))
						html_w.write('''</ul></li>''')
					if no_num > 0:
						html_w.write('''<li style=\\"color: #CD6600;\\">Protein no regulated <ul>''')
						for KO in KOs.split('+'):
							for pro in KO2pro[KO].split(';'):
								if pro in pro2color.keys():
									if 'yellow' == pro2color[pro]:
										FC = pro2FC[pro]
										html_w.write('''<li>%s: %s (%s) </li>''' %(KO,pro,FC))
						html_w.write('''</ul></li>''')
					
					if up_num >0 and down_num == 0:
						'''只有上调和无明显表达'''
						x1_p = x1 + (x2 - x1)/2
						y1_p = y1 + (y2 - y1)/2
						x2_p = x2
						y2_p = y2
						x3_p = x3
						y3_p = y3
						x4_p = x3 + (x4 - x3)/2
						y4_p = y3 + (y4 - y3)/2
						
						layer.polygon([(x1_p,y1_p),(x2_p,y2_p),(x3_p,y3_p),(x4_p,y4_p)], fill=(255, 0, 0, 128))
						
					if up_num == 0 and down_num > 0:
						'''只有下调和无明显表达'''
						x1_p = x1 + (x2 - x1)/2
						y1_p = y1 + (y2 - y1)/2
						x2_p = x2
						y2_p = y2
						x3_p = x3
						y3_p = y3
						x4_p = x3 + (x4 - x3)/2
						y4_p = y3 + (y4 - y3)/2
						
						layer.polygon([(x1_p,y1_p),(x2_p,y2_p),(x3_p,y3_p),(x4_p,y4_p)], fill=(0, 255, 0, 128))
						
					if up_num == 0 and down_num == 0 and no_num > 0:
						'''只有无明显表达'''
						x1_p = x1 + (x2 - x1)/2
						y1_p = y1 + (y2 - y1)/2
						x2_p = x2
						y2_p = y2
						x3_p = x3
						y3_p = y3
						x4_p = x3 + (x4 - x3)/2
						y4_p = y3 + (y4 - y3)/2
						
						layer.polygon([(x1_p,y1_p),(x2_p,y2_p),(x3_p,y3_p),(x4_p,y4_p)], fill=(255, 255, 0, 128))
						
					if up_num > 0 and down_num > 0 and no_num == 0:
						'''只有上调和下调'''
						x1_p = x1 + (x2 - x1)/2
						y1_p = y1 + (y2 - y1)/2
						x2_p = x2
						y2_p = y2
						x5_p = x3
						y5_p = y3
						x6_p = x3 + (x4 - x3)/2
						y6_p = y3 + (y4 - y3)/2
						x3_p = x1_p + (x5_p - x1_p)/2
						y3_p = y1_p + (y5_p - y1_p)/2
						x4_p = x2_p + (x6_p - x2_p)/2
						y4_p = y2_p + (y6_p - y2_p)/2
						
						layer.polygon([(x1_p,y1_p),(x2_p,y2_p),(x4_p,y4_p),(x3_p,y3_p)], fill=(255, 0, 0, 128))
						layer.polygon([(x3_p,y3_p),(x4_p,y4_p),(x6_p,y6_p),(x5_p,y5_p)], fill=(0, 255, 0, 128))						
						
					if up_num > 0 and down_num > 0 and no_num > 0:
						'''上调下调无明显表达都有'''		
						x1_p = x1 + (x2 - x1)/2
						y1_p = y1 + (y2 - y1)/2
						x2_p = x2
						y2_p = y2
						x7_p = x3 + (x4 - x3)/2
						y7_p = y3 + (y4 - y3)/2
						x8_p = x3
						y8_p = y3
						x3_p = x1_p + (x7_p - x1_p)/3*1
						y3_p = y1_p + (y7_p - y1_p)/3*1
						x4_p = x2_p + (x8_p - x2_p)/3*1
						y4_p = y2_p + (y8_p - y2_p)/3*1
						x5_p = x1_p + (x7_p - x1_p)/3*2
						y5_p = y1_p + (y7_p - y1_p)/3*2
						x6_p = x2_p + (x8_p - x2_p)/3*2
						y6_p = y2_p + (y8_p - y2_p)/3*2
						
						layer.polygon([(x1_p,y1_p),(x2_p,y2_p),(x4_p,y4_p),(x3_p,y3_p)], fill=(255, 0, 0, 128))
						layer.polygon([(x3_p,y3_p),(x4_p,y4_p),(x6_p,y6_p),(x5_p,y5_p)], fill=(255, 255, 0, 128))
						layer.polygon([(x5_p,y5_p),(x6_p,y6_p),(x8_p,y8_p),(x7_p,y7_p)], fill=(0, 255, 0, 128))
				
			final = Image.new("RGBA", background.size)
			final.paste(background, (0,0), background)
			final.paste(coord_draw, (0,0), coord_draw)
			final.save('%s/%s.png' %(outdir,name))
			html_w.write('''</ul>");' />
''')	
	
	for KO in set(KOs_set.split('+')):
		if KO in KO2rna.keys() or KO in KO2pro.keys():
			path_combine += '%s\t%s\t' %(name,KO)
			if KO in KO2rna.keys():
				for rna in KO2rna[KO].split(';'):
					if rna in rna2FC.keys():
						path_combine += '(%s;%s);' %(rna,rna2FC[rna])
			path_combine += '\t'
			if KO in KO2pro.keys():
				for pro in KO2pro[KO].split(';'):
					if pro in pro2FC.keys():
						path_combine += '(%s;%s);' %(pro,pro2FC[pro])
			path_combine += '\n'
	
	html_w.write('''
</map>
<img src='./%s.png' usemap='#%s.png' />
<div id='result' style='position: absolute; width: 50%%; border: 1px solid #000; background-color: #fff; filter: alpha(opacity=95); opacity: 0.95; font-size: 12px; padding-right: 20px; display: none;' onmouseover="javascript: this.style.filter = 'alpha(opacity=100)'; this.style.opacity = 1;" onmouseout="javascript: this.style.filter ='alpha(opacity=95)'; this.style.opacity = 0.95;"></div>
</body></html>'''%(name,name))	
	html_w.close()
	return path_combine

with open ('%s/Combined-KO-path.xls' %args.out,'w')as path_w:
	path_w.write('Paths\tKO\tmRNAs\tProteins\n')	
	for path in os.listdir(args.path):
		if 'png' in path:
			name = path.split('.png')[0]		
			com_xls = MRNA_pro_KEGG(name,args.path,args.out,KO2rna,rna2color,rna2FC,KO2pro,pro2color,pro2FC)
			path_w.write(str(com_xls))