#!/usr/bin/env python

from __future__ import division
from datetime import datetime
from subprocess import Popen, PIPE
from optparse import OptionParser
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import scipy.stats as st
import copy,os,glob,sys,random,math,warnings,ast,time,pickle,shutil,operator,scipy,re,cProfile
matplotlib.rcParams['figure.figsize'] = (14, 10.5)
matplotlib.rcParams['font.size'] = 18
warnings.simplefilter('ignore')
# matplotlib.style.use('ggplot')


def make_directory(folder,save_data):
	output_folder=(folder+'/')
	x=1
	while True:
		if not os.path.isdir(output_folder): 
			try:
				os.makedirs(output_folder)
			except Exception as e:
				pass
			else:
				os.chdir(output_folder)
				os.makedirs("temp")
				if save_data:
					os.makedirs("history_data")
				break
		else:
			output_folder=(folder+'-'+str(x)+'/')
			if not os.path.isdir(output_folder):
				try:
					os.makedirs(output_folder)
				except Exception as e:
					pass
				else:
					os.chdir(output_folder)
					os.makedirs("temp")
					if save_data:
						os.makedirs("history_data")
					break
			else:
				x+=1
		


def data_process(bam1,bam2,lib,case):
	samtool="/share/apps/samtools/bin/samtools"   # Please put your samtools path here
	try:
		x=Popen([samtool,"--help"], stdout=PIPE, stderr=PIPE).communicate()
	except Exception as e:
		sys.exit("Please put samtools path at line 53")
	
	def split_chromosome(all_chromosome,bam1,bam2):
		for chromosome in all_chromosome:
			count_read['1']["+"].update({chromosome:{}})
			count_read['1']["-"].update({chromosome:{}})
			count_read['2']["+"].update({chromosome:{}})
			count_read['2']["-"].update({chromosome:{}})
			print_process("calculate depth "+lib+" : chromosome : "+chromosome+"...",True)
			data, err=Popen([samtool,"view","../"+bam1,chromosome], stdout=PIPE, stderr=PIPE).communicate()
			data=data.split('\n')
			data=data[:-1]
			count_read['1']["+"][chromosome],count_read['1']["-"][chromosome],total_reads=find_depth(
									data,chromosome,case)
			count_read['total_reads_1']+=total_reads
			data, err=Popen([samtool,"view","../"+bam2,chromosome], stdout=PIPE, stderr=PIPE).communicate()
			data=data.split('\n')
			data=data[:-1]
			count_read['2']["+"][chromosome],count_read['2']["-"][chromosome],total_reads=find_depth(
									data,chromosome,case)
			count_read['total_reads_2']+=total_reads
			if memory:
				save_varible([count_read['1']["+"][chromosome],count_read['2']["+"][chromosome]],"temp/"+lib+"/+/depth_"+chromosome)
				count_read['1']["+"][chromosome]={}
				count_read['2']["+"][chromosome]={}
				save_varible([count_read['1']["-"][chromosome],count_read['2']["-"][chromosome]],"temp/"+lib+"/-/depth_"+chromosome)
				count_read['1']["-"][chromosome]={}
				count_read['2']["-"][chromosome]={}		
		return count_read
	header,err=Popen([samtool,"view","-H","../"+bam1], stdout=PIPE, stderr=PIPE).communicate()
	header2,err=Popen([samtool,"view","-H","../"+bam2], stdout=PIPE, stderr=PIPE).communicate()
	if err != "" :
		sys.exit(err)
	elif header == "":
		sys.exit("the header of "+bam1+" missing")
	elif header2 == "":
		sys.exit("the header of "+bam2+" missing")
	header=header.split('\n')
	all_chromosome={}
	for line in header:
		if '@HD' in line:
			issort=line.split("\t")[-1].replace("SO:","")
		elif '@SQ' in line:
			line=line.split('\t')
			chromosome = line[1].replace('SN:','')
			length =int(line[2].replace('LN:',''))
			all_chromosome.update({chromosome:length})
	header2=header2.split('\n')
	for line in header2:
		if '@HD' in line:
			issort2=line.split("\t")[-1].replace("SO:","")
			break
	if issort == "coordinate":
		if not os.path.isfile("../"+bam1+".bai"):
			print_process("create index of " + bam1 +" file...",True)
			not_use=Popen([samtool,"index","../"+bam1], stdout=PIPE, stderr=PIPE).communicate()
	else:
		print_process("sort "+bam1+" file...",True)
		not_use=Popen([samtool,"sort","-T","/tmp/"+bam1.replace(".bam","_sort.sorted"),"-o","../"+bam1.replace(".bam","_sort.bam"),"../"+bam1], stdout=PIPE, stderr=PIPE).communicate()
		os.remove("../"+bam1)
		os.rename("../"+bam1.replace(".bam","_sort.bam"),"../"+bam1)
		print_process("create index of " + bam1 +" file...",True)
		not_use=Popen([samtool,"index","../"+bam1], stdout=PIPE, stderr=PIPE).communicate()
	if issort2 == "coordinate":
		if not os.path.isfile("../"+bam2+".bai"):
			print_process("create index of " + bam2 +" file...",True)
			not_use=Popen([samtool,"index","../"+bam2], stdout=PIPE, stderr=PIPE).communicate()
	else:
		print_process("sort "+bam1+" file...",True)
		not_use=Popen([samtool,"sort","-T","/tmp/"+bam2.replace(".bam","_sort.sorted"),"-o","../"+bam2.replace(".bam","_sort.bam"),"../"+bam2], stdout=PIPE, stderr=PIPE).communicate()
		os.remove("../"+bam2)
		os.rename("../"+bam2.replace(".bam","_sort.bam"),"../"+bam2)
		print_process("create index of " + bam2 +" file...",True)
		not_use=Popen([samtool,"index","../"+bam2], stdout=PIPE, stderr=PIPE).communicate()
	count_read={'1':{"+":{},"-":{}},'2':{"+":{},"-":{}},'total_reads_1':0,'total_reads_2':0}
	count_read = split_chromosome(all_chromosome,bam1,bam2)
	return count_read,all_chromosome

def find_depth(input_reads,chromosome,case):
	position_reads={'+':{},'-':{}}
	position_reads['+'].update({chromosome:{}})
	position_reads['-'].update({chromosome:{}})
	num_position=0
	for line in input_reads:
		inline=line.split('\t')
		try:
			start=int(inline[3])
			chromosome=inline[2]
			CIGAR=inline[5]
			seq=inline[9]
			len_seq=len(seq)
			read_quality=int(inline[4])
			base_score=inline[10]
		except Exception as e:
			pass
		try:
			FLAG=bin(int(inline[1]))
			FLAG=FLAG[2:]
			if len(FLAG) < 12:
				for x in range(0,12-len(FLAG)):
					FLAG='0'+FLAG
		except:
			FLAG='000000000000'
		if FLAG[-5] == '0' and FLAG[-3] != "1": #strand forward and not unmapped read
			last=0
			num_seq=0
			for x in range(0,len(CIGAR)):
				if not CIGAR[x].isdigit():	
					if case == 'start':	
						if CIGAR[x] == 'M': #alignment match
							try:
								position_reads['+'][chromosome][start] +=1
								num_position+=1
							except:
								position_reads['+'][chromosome].update({start:1})
								num_position+=1
							break	
						elif CIGAR[x] == 'I' or CIGAR[x] == 'S': #Insertion to reference (no position in reference)
							num_seq+=int(CIGAR[last:x])
							last=x+1
						elif CIGAR[x] == 'D' or CIGAR[x] == 'N' :	
							try:
								start=start+int(CIGAR[last:x])
							except:
								pass
							last=x+1
						else:
							last=x+1
					elif case == "all":	#whole read
						if CIGAR[x] == 'M':#alignment match
							for position in range(start,start+int(CIGAR[last:x])):
								try:
									position_reads['+'][chromosome][position] +=1
									num_position+=1
								except:
									position_reads['+'][chromosome].update({position:1})
									num_position+=1
								num_seq+=1
							start=start+int(CIGAR[last:x])
							last=x+1							
						elif CIGAR[x] == 'I' or CIGAR[x] == 'S': #Insertion to reference (no position in reference)
							num_seq+=int(CIGAR[last:x])
							last=x+1
						elif CIGAR[x] == 'D' or CIGAR[x] == 'N' :	
							try:
								start=start+int(CIGAR[last:x])
							except:
								pass
							last=x+1
						else:
							last=x+1
					elif case == "end":
						if CIGAR[x] == 'M' :
							start=start+int(CIGAR[last:x])
							last=x+1
						elif CIGAR[x] == 'I' or CIGAR[x] == 'S': #Insertion to reference (no position in reference)
							num_seq+=int(CIGAR[last:x])
							last=x+1
						elif CIGAR[x] == 'D' or CIGAR[x] == 'N' :	
							try:
								start=start+int(CIGAR[last:x])
							except:
								pass
							last=x+1
						else:
							last=x+1
			if case == 'end':
				try:
					position_reads['+'][chromosome][start-1] +=1
					num_position+=1
				except:
					position_reads['+'][chromosome].update({start-1:1})
					num_position+=1		
		elif FLAG[-5] == '1' and FLAG[-3] != "1": 	#strand reverse and not unmapped read
			last=0
			num_seq=0
			for x in range(0,len(CIGAR)):
				if not CIGAR[x].isdigit():
					if case == 'start': ## Start_reads
						if CIGAR[x] == 'M' :
							start=start+int(CIGAR[last:x])
							last=x+1
						elif CIGAR[x] == 'I' or CIGAR[x] == 'S': #Insertion to reference (no position in reference)
							num_seq+=int(CIGAR[last:x])
							last=x+1
						elif CIGAR[x] == 'D' or CIGAR[x] == 'N' :	
							try:
								start=start+int(CIGAR[last:x])
							except:
								pass
							last=x+1	
						else:
							last=x+1
					elif case == 'end':	
						if CIGAR[x] == 'M': #alignment match
							try:
								position_reads['-'][chromosome][start] +=1
								num_position+=1
							except:
								position_reads['-'][chromosome].update({start:1})
								num_position+=1
							break	
						elif CIGAR[x] == 'I' or CIGAR[x] == 'S': #Insertion to reference (no position in reference)
							num_seq+=int(CIGAR[last:x])
							last=x+1
						elif CIGAR[x] == 'D' or CIGAR[x] == 'N' :	
							try:
								start=start+int(CIGAR[last:x])
							except:
								pass
							last=x+1
						else:
							last=x+1
					elif case == 'all': ## All-reads
						if CIGAR[x] == 'M':
							for position in range(start,start+int(CIGAR[last:x])):
								try:
									position_reads['-'][chromosome][position] +=1
									num_position+=1
								except:
									position_reads['-'][chromosome].update({position:1})
									num_position+=1	
								num_seq+=1
							start=start+int(CIGAR[last:x])
							last=x+1
						elif CIGAR[x] == 'I' or CIGAR[x] == 'S': #Insertion to reference (no position in reference)
							num_seq+=int(CIGAR[last:x])
							last=x+1
						elif CIGAR[x] == 'D' or CIGAR[x] == 'N' :	
							try:
								start=start+int(CIGAR[last:x])
							except:
								pass
							last=x+1
						else:
							last=x+1
			if case == 'start':
				try:
					position_reads['-'][chromosome][start-1] +=1
					num_position+=1
				except:
					position_reads['-'][chromosome].update({start-1:1})
					num_position+=1
	return [position_reads['+'][chromosome],position_reads['-'][chromosome],num_position]

def check_pseudo(input_depth_1,input_depth_2,reads_fillter,filter_type):
	ratio_all_position={}
	both=0
	only_1=0
	only_2=0
	for strand in input_depth_1:
		ratio_all_position.update({strand:{}})
		for chromosome in sorted(input_depth_1[strand]):
			ratio_all_position[strand].update({chromosome:{}})
			for position in sorted(input_depth_1[strand][chromosome]):
				depth_1=input_depth_1[strand][chromosome][position]
				try:
					depth_2=input_depth_2[strand][chromosome][position]
				except Exception as e:
					if depth_1 >= reads_fillter:
						ratio_all_position[strand][chromosome].update({position:True})
						only_1+=1
				else:
					if filter_type == 'total':
						if depth_1+depth_2 >= reads_fillter:
							ratio_all_position[strand][chromosome].update({position:True})
							both+=1
					else:
						if depth_1 >= reads_fillter and depth_2 >= reads_fillter:
							ratio_all_position[strand][chromosome].update({position:True})
							both+=1					
	for strand in input_depth_2:
		for chromosome in sorted(input_depth_2[strand]):
			x=set(input_depth_2[strand][chromosome].keys())-set(ratio_all_position[strand][chromosome].keys())
			for position in x:
				try:
					depth_1=input_depth_1[strand][chromosome][position]
				except:
					if input_depth_2[strand][chromosome][position] >= reads_fillter:
						depth_2=input_depth_2[strand][chromosome][position]
						only_2+=1
	return [both,only_1,only_2]
def find_ratio(input_depth_1,num_1,input_depth_2,num_2,pseudo,lib,save_data,reads_count,filter_type):
	if save_data:
		write_depth1=open('history_data/'+lib+"_depth_1.txt",'w')
		write_depth1.write('\t'.join(['strand','chromosome','position','depth'])+'\n')
		write_depth2=open('history_data/'+lib+"_depth_2.txt",'w')
		write_depth2.write('\t'.join(['strand','chromosome','position','depth'])+'\n')
	ratio_all_position={}
	both=0
	only_1=0
	only_2=0
	num_normalize=float(num_2/num_1)
	# data1=[]
	# data2=[]
	for strand in sorted(input_depth_1):
		ratio_all_position.update({strand:{}})
		for chromosome in sorted(input_depth_1[strand]):
			ratio_all_position[strand].update({chromosome:{}})
			if memory:
				input_depth_1[strand][chromosome],input_depth_2[strand][chromosome]=load_varible("temp/"+lib+"/"+strand+"/depth_"+chromosome)
			for position in sorted(input_depth_1[strand][chromosome]):
				depth_1=input_depth_1[strand][chromosome][position]
				# data1.append(depth_1)
				if not pseudo:
					try:
						depth_2=input_depth_2[strand][chromosome][position]
					except Exception as e:
						if depth_1 >= reads_count:
							only_1+=1
							if save_data:
								write_depth1.write('\t'.join([strand,chromosome,str(position),str(depth_1)])+'\n')
					else:
						# data2.append(depth_2)
						if filter_type == 'total':
							if depth_1+depth_2 >= reads_count:
								ratio_at_position=float(depth_1/depth_2)*num_normalize
								ratio_all_position[strand][chromosome].update({position:ratio_at_position})
								both+=1
								if save_data:
									write_depth1.write('\t'.join([strand,chromosome,str(position),str(depth_1)])+'\n')
									write_depth2.write('\t'.join([strand,chromosome,str(position),str(depth_2)])+'\n')
						else:
							if depth_1 >= reads_count and depth_2 >= reads_count:
								ratio_at_position=float(depth_1/depth_2)*num_normalize
								ratio_all_position[strand][chromosome].update({position:ratio_at_position})
								both+=1
								if save_data:
									write_depth1.write('\t'.join([strand,chromosome,str(position),str(depth_1)])+'\n')
									write_depth2.write('\t'.join([strand,chromosome,str(position),str(depth_2)])+'\n')
				else:
					try:
						depth_2=input_depth_2[strand][chromosome][position]
					except Exception as e:
						if depth_1 >= reads_count:
							ratio_at_position=float(depth_1+1)*num_normalize
							ratio_all_position[strand][chromosome].update({position:ratio_at_position})
							only_1+=1
							if save_data:
								write_depth1.write('\t'.join([strand,chromosome,str(position),str(depth_1)])+'\n')
					else:
						# data2.append(depth_2)
						if filter_type == 'total':
							if depth_1+depth_2 >= reads_count:
								ratio_at_position=float((depth_1+1)/(depth_2+1))*num_normalize
								ratio_all_position[strand][chromosome].update({position:ratio_at_position})
								both+=1
								if save_data:
									write_depth1.write('\t'.join([strand,chromosome,str(position),str(depth_1)])+'\n')
									write_depth2.write('\t'.join([strand,chromosome,str(position),str(depth_2)])+'\n')
						else:
							if depth_1 >= reads_count and depth_2 >= reads_count:
								ratio_at_position=float((depth_1+1)/(depth_2+1))*num_normalize
								ratio_all_position[strand][chromosome].update({position:ratio_at_position})
								both+=1
								if save_data:
									write_depth1.write('\t'.join([strand,chromosome,str(position),str(depth_1)])+'\n')
									write_depth2.write('\t'.join([strand,chromosome,str(position),str(depth_2)])+'\n')
			if memory:
				save_varible([input_depth_1[strand][chromosome],input_depth_2[strand][chromosome]],"temp/"+lib+"/"+strand+"/depth_"+chromosome)
				save_varible(ratio_all_position[strand][chromosome],"temp/"+lib+"/"+strand+"/ratio_"+chromosome)
				input_depth_1[strand][chromosome]={}
				input_depth_2[strand][chromosome]={}
				ratio_all_position[strand][chromosome]={}
	for strand in input_depth_2:
		for chromosome in sorted(input_depth_2[strand]):
			if memory:
				input_depth_1[strand][chromosome],input_depth_2[strand][chromosome]=load_varible("temp/"+lib+"/"+strand+"/depth_"+chromosome)
				ratio_all_position[strand][chromosome]=load_varible("temp/"+lib+"/"+strand+"/ratio_"+chromosome)
			x=set(input_depth_2[strand][chromosome].keys())-set(ratio_all_position[strand][chromosome].keys())
			for position in x:
				depth_2=input_depth_2[strand][chromosome][position]
				try:
					depth_1=input_depth_1[strand][chromosome][position]
				except:
					# data2.append(depth_2)
					if depth_2 >= reads_count:
						only_2+=1
						if save_data:
							write_depth2.write('\t'.join([strand,chromosome,str(position),str(depth_2)])+'\n')
						if pseudo:
							ratio_at_position=float(1/(depth_2+1))*num_normalize
							ratio_all_position[strand][chromosome].update({position:ratio_at_position})
			if memory:
				save_varible([input_depth_1[strand][chromosome],input_depth_2[strand][chromosome]],"temp/"+lib+"/"+strand+"/depth_"+chromosome)
				save_varible(ratio_all_position[strand][chromosome],"temp/"+lib+"/"+strand+"/ratio_"+chromosome)
				input_depth_1[strand][chromosome]={}
				input_depth_2[strand][chromosome]={}
				ratio_all_position[strand][chromosome]={}
	if save_data:
		write_depth1.close()
		write_depth2.close()
	return [ratio_all_position,both,only_1,only_2]	
	# return [ratio_all_position,both,only_1,only_2,data1,data2]

def boxcox_transform(enrichment,lamda):
	if lamda != 0:
		enrichment_tranform = (enrichment**lamda - 1) / lamda
	else:
		enrichment_tranform = np.log(enrichment)
	return enrichment_tranform
def boxcox_reverse(enrichment,lamda):
	enrichment_reverse=np.exp(np.log(lamda*enrichment+1)/lamda)
	return enrichment_reverse
def boxcox_data_transform(ratio,lib,p_value,pass_cut_off):
	data=[]
	for strand in ratio:
		for chromosome in ratio[strand]:
			if memory:
				ratio[strand][chromosome]=load_varible("temp/"+lib+"/"+strand+"/ratio_"+chromosome)
			data+=ratio[strand][chromosome].values()
			if memory:
				ratio[strand][chromosome]={}
	print_out("\tall_position = "+str(len(data)),True)
	takeClosest = lambda num,collection:min(collection,key=lambda x:abs(x-num))
	sp=plt.subplot(2, 1,1)
	pp=st.probplot(data, dist=st.norm, plot=sp)
	axes = plt.gca()
	y_min,y_max=axes.get_ylim()
	x_min,x_max=axes.get_xlim()
	r_before=pp[1][2]
	plt.text(0,y_max*0.8,r'$R^2$'+' = '+str(round(r_before**2,4)),fontsize=18)
	plt.title('Before Box-Cox transformation',fontsize=18)
	sp=plt.subplot(2, 1,2)
	boxcox,lamda=st.boxcox(data)
	pp=st.probplot(boxcox, dist=st.norm, plot=sp)
	y_p=pp[0][1].tolist()
	r_after=pp[1][2]
	params = st.norm.fit(boxcox)
	arg = params[:-2]
	loc = params[-2]
	scale = params[-1]
	y=st.norm.ppf(1-p_value, loc=loc, scale=scale, *arg)
	if np.isfinite(y):
		i=y_p.index(takeClosest(y,y_p))
	elif np.isinf(y) and y < 0:
		y=min(boxcox)
		i=1
	elif np.isinf(y) and y > 0:
		y=max(boxcox)
		i=len(boxcox)
	plt.axhline(y=y,color='g',label="p_value "+str(p_value))
	axes = plt.gca()
	y_min,y_max=axes.get_ylim()
	x_min,x_max=axes.get_xlim()
	center_dot=y_p[int(len(y_p)/2)]
	center_graph=(y_max-y_min)/2
	y_old=boxcox_reverse(y,lamda)
	if center_dot >= center_graph:
		height=(y_max-center_dot)/4
		plt.text(-2,center_graph+(3*height),'lambda = '+str(round(lamda,4)),fontsize=18)
		plt.text(-2,center_graph+(height),"cutoff value = "+str(round(y_old,2))+"("+str(round(y,2))+")",fontsize=18)
	else:
		height=(center_dot-y_min)/4
		plt.text(-2,center_graph-(height),'lambda = '+str(round(lamda,4)),fontsize=18)
		plt.text(-2,center_graph-(3*height),"cutoff value = "+str(round(y_old,2))+"("+str(round(y,2))+")",fontsize=18)
	plt.tight_layout()
	plt.legend(loc='upper left', shadow=True,fontsize=16)
	plt.title('After Box-Cox transformation',fontsize=18)
	plt.savefig(lib+'_Box-Cox_transformation.png', dpi=300)
	plt.clf()
	if r_after**2 > pass_cut_off:
		return [lamda,params,y_old,y,r_after**2,True]
	else:
		return [lamda,params,y_old,y,r_after**2,False]

def save_history_data(ratio,all_lamda,lib,save_data,boxcox):
	if boxcox:
		lamda=all_lamda[lib]
	if save_data:
		write_ratio=open('history_data/'+lib+"_ratio.txt",'w')
		if boxcox:
			write_ratio.write('\t'.join(['strand','chromosome','position','ratio(Before boxcox)','ratio(After boxcox)'])+'\n')
		else:
			write_ratio.write('\t'.join(['strand','chromosome','position','ratio'])+'\n')			
	have_to={}
	for strand in sorted(ratio):
		have_to.update({strand:{}})
		for chromosome in sorted(ratio[strand]):
			have_to[strand].update({chromosome:{}})
			if memory:
				ratio[strand][chromosome]=load_varible("temp/"+lib+"/"+strand+"/ratio_"+chromosome)
			for position in sorted(ratio[strand][chromosome]):
				enrichment=ratio[strand][chromosome][position]
				if boxcox:
					try:
						enrichment_tranform=have_to[strand][chromosome][position]
					except:
						enrichment_tranform=boxcox_transform(enrichment,lamda)
						have_to[strand][chromosome].update({position:enrichment_tranform})
					# ratio[strand][chromosome][position]=enrichment_tranform
					if save_data:
						write_ratio.write('\t'.join([strand,chromosome,str(position),str(enrichment),str(enrichment_tranform)])+'\n')
				else:
					if save_data:
						write_ratio.write('\t'.join([strand,chromosome,str(position),str(enrichment)])+'\n')
			if memory:
				# save_varible(ratio[strand][chromosome],"temp/"+lib+"/"+strand+"/ratio_"+chromosome)
				ratio[strand][chromosome]={}
	if save_data:
		write_ratio.close()

def find_empirical(ratio,percentile,lib):
	takeClosest = lambda num,collection:min(collection,key=lambda x:abs(x-num))
	data=[]
	for strand in ratio:
		for chromosome in ratio[strand]:
			if memory:
				ratio[strand][chromosome]=load_varible("temp/"+lib+"/"+strand+"/ratio_"+chromosome)
			data+=ratio[strand][chromosome].values()
			if memory:
				ratio[strand][chromosome]={}
	print_out("\tall_position = "+str(len(data)),False)
	data_unique=list(set(data))
	data_unique.sort()
	density = st.gaussian_kde(data)
	xs = np.linspace(min(data),max(data),100)
	density.covariance_factor = lambda : .25
	density._compute_covariance()
	plt.plot(xs,density(xs),label="empirical distribution")
	plt.fill_between(xs,0, density(xs), facecolor='blue', alpha=0.5)
	peak=st.scoreatpercentile(data,percentile)
	peak_index=data_unique.index(takeClosest(peak,data_unique))
	if st.percentileofscore(data,peak) < percentile:
		while True:
			peak_index=data_unique.index(takeClosest(peak,data_unique))
			peak=data_unique[peak_index+1]
			if st.percentileofscore(data,peak) >= percentile: break
	axes = plt.gca()
	y_min,y_max=axes.get_ylim()
	x_min,x_max=axes.get_xlim()
	plt.axvline(x=peak,color='r')
	plt.text(peak,y_max*0.6,"percentile = "+str(percentile)+", critical value = "+str(round(peak,5)),fontsize=18)
	plt.title(lib)
	plt.xlabel('Ratio')
	plt.ylabel('Probability density')
	plt.legend(loc='upper right', shadow=True)
	plt.savefig(lib+'_top-rank.png', dpi=300)
	plt.clf()
	all_p_value={}	
	data_unique=data_unique[peak_index:]
	return peak


def pie_plot(data,text,name_save,lib):
	def make_autopct(values):
		def my_autopct(pct):
			total = sum(values)
			val = int(round(pct*total/100.0))
			return '{p:.2f}%  ({v:d})'.format(p=pct,v=val)
		return my_autopct
	labels = text
	sizes = data
	color = ['salmon','deepskyblue','lightgreen']
	explode=[0]*len(data)
	explode[0]=0.2
	explode = tuple(explode)  # explode 1st slice
	plt.pie(sizes, labels=labels, colors=color,autopct=make_autopct(sizes),explode=explode, startangle=140)
	plt.axis('equal')
	plt.title(lib)
	plt.savefig(name_save+'.png', dpi=300)
	plt.clf()


def input_gene_gff(input_folder):
	takeClosest = lambda num,collection:min(collection,key=lambda x:abs(x-num))
	gene_format={'+':{},'-':{}}
	gff_file=glob.glob(os.path.join('../', '*.gff'))
	if len(gff_file) == 0: 
		print_out("\ngene annotate :\n\tcouldn't put gene annotate to output file because No such gff file in "+input_folder,False)
		return False,False
	elif len(gff_file) > 1:
		print_out("\ngene annotate :\n\tcouldn't put gene annotate to output file because there are more than one gff file in "+input_folder,False)
		return False,False
	f=open(gff_file[0])
	f=f.readlines()
	have_chromosome=False
	new_chromosome=''
	region_now={'+':{},'-':{}}
	gene_region={'+':{},'-':{}}
	for line in f:
		line=line.replace('\n','')
		if line[0]=='#':
			if '##FASTA' in line: break	
			if '##sequence-region' in line:
				have_chromosome=True
				inline=line.split()
				len_chromosome=int(inline[3])
				chromosome=inline[1]
				gene_format['+'].update({chromosome:{'length':len_chromosome,'annotate':{}}})
				gene_format['-'].update({chromosome:{'length':len_chromosome,'annotate':{}}})
				region_now['+'].update({chromosome:[0,0]})
				region_now['-'].update({chromosome:[0,0]})
				gene_region['+'].update({chromosome:[[]]})
				gene_region['-'].update({chromosome:[[]]})
		else:
			if not have_chromosome:
				print_out("\ngene annotate :\n\tcouldn't put gene annotate to output file because gff file don't have header line '##sequence-region'",False)
				return False,False
			inline=line.split('\t')
			feature=inline[2]
			chromosome=inline[0]
			start=int(inline[3])
			end=int(inline[4])
			size=(abs(end-start))
			strand=inline[6]
			info=inline[-1].split(';')
			gene_name=""
			locus_tag=""
			for text in info:
				if re.search('name=',text,re.IGNORECASE):
					try:
						gene_name=text.split('=')[-1]
					except:
						pass
				elif re.search('locus',text,re.IGNORECASE):
					try:
						locus_tag=text.split('=')[-1]
					except:
						pass
			if not gene_name:
				if locus_tag:
					gene_name=locus_tag
				else:
					gene_name=str(start)
			if not locus_tag:
				locus_tag="n/a"
			if start !=1 or (start ==1 and end != gene_format[strand][chromosome]['length']):
				all_position=sorted(gene_format[strand][chromosome]['annotate'])
				try:
					near=takeClosest(start,all_position)
				except:
					gene_format[strand][chromosome]['annotate'].update({start:{'name':gene_name,'end':end,'locus':locus_tag,'inside':{},'feature':feature}})
				else:
					index_near=all_position.index(near)
					position_detail=gene_format[strand][chromosome]['annotate'][near]
					if start==near and end == position_detail['end']:
						gene_format[strand][chromosome]['annotate'][near]['feature']+=','+feature
						gene_format[strand][chromosome]['annotate'][near]['name']+=','+gene_name
					elif near <= start <= position_detail['end'] and near <= end <= position_detail['end']:
						if start in position_detail['inside']:
							gene_format[strand][chromosome]['annotate'][near]['inside'][start]['name']+=','+gene_name
							gene_format[strand][chromosome]['annotate'][near]['inside'][start]['feature']+=','+feature
						else:	
							gene_format[strand][chromosome]['annotate'][near]['inside'].update({start:{'end':end,'name':gene_name,'feature':feature}})
					elif start <= near and end >= position_detail['end']:
						
						gene_format[strand][chromosome]['annotate'].update({start:{'name':gene_name,'end':end,'locus':locus_tag,'feature':feature,'inside':{
								near:{'end':position_detail['end'],'name':position_detail['name'],'feature':position_detail['feature']}
							}}})
						if start!=near:
							del gene_format[strand][chromosome]['annotate'][near]
					elif (start < near and end < near) or (start > position_detail['end'] and end > position_detail['end']):
						if index_near == 0 or index_near == len(all_position)-1:
							gene_format[strand][chromosome]['annotate'].update({start:{'name':gene_name,'end':end,'locus':locus_tag,'inside':{},'feature':feature}})
						else:
							if start > near:
								near=all_position[index_near+1]
							else:
								near=all_position[index_near-1]	
							position_detail=gene_format[strand][chromosome]['annotate'][near]
							if near <= start <= position_detail['end'] and near <= end <= position_detail['end']:
								if start in position_detail['inside']:
									gene_format[strand][chromosome]['annotate'][near]['inside'][start]['name']+=','+gene_name
									gene_format[strand][chromosome]['annotate'][near]['inside'][start]['feature']+=','+feature
								else:	
									gene_format[strand][chromosome]['annotate'][near]['inside'].update({start:{'end':end,'name':gene_name,'feature':feature}})
							else:
								gene_format[strand][chromosome]['annotate'].update({start:{'name':gene_name,'end':end,'locus':locus_tag,'inside':{},'feature':feature}})
					else:
						# print strand,chromosome,feature,start,end,near,position_detail['end']
						gene_format[strand][chromosome]['annotate'].update({start:{'name':gene_name,'end':end,'locus':locus_tag,'inside':{},'feature':feature}})
	num_gene=0
	for strand in gene_format:
		for chromosome in gene_format[strand]:
			num_gene+=len(gene_format[strand][chromosome]['annotate'])
	return [gene_format,num_gene]

def write_defult_result(ratio,depth,p_value_cut_off,params_all,peak_cut_off_all,gene_format,annotate,boxcox_lamda,distribution):	
	strand_check=['+','-']
	all_gene_position_sort={}
	if gene_format:
		for strand in gene_format:
			all_gene_position_sort.update({strand:{}})
			for chromosome in gene_format[strand]:
				all_gene_position_sort[strand].update({chromosome:gene_format[strand][chromosome]['annotate'].keys()})
				all_gene_position_sort[strand][chromosome].sort()
	for lib in sorted(ratio):
		if distribution == 'toprank':
			lib_sig=True
		else:
			if lib in boxcox_lamda:
				lib_sig=True
			else:
				lib_sig=False
		write_file=open(lib+'_result.txt','w')
		if lib_sig:
			peak_cut_off=peak_cut_off_all[lib]
			print_process("write result "+lib+"...",True)
			write_gff={}
			if gene_format:
				write_file.write('position\tstrand\tchromosome\tdistance\tregion\tfeature\tf_name\tf_strand\tf_start\tf_end\t')
			else:
				write_file.write('position\tstrand\tchromosome\t')
			if distribution == 'toprank':
				write_file.write('\t'.join(['depth_1','depth_2','ratio','\n']))
			else:
				lamda=boxcox_lamda[lib]
				write_file.write('\t'.join(['depth_1','depth_2','ratio(before Box-Cox)','ratio(after Box-Cox)','p-value','\n']))
			write_gff={'+':open(lib+'_forward.gff','w'),'-':open(lib+'_reverse.gff','w')}	
			write_gff['+'].write("#gff-version 3\n")	
			write_gff['-'].write("#gff-version 3\n")			
			for strand in ratio[lib]:
				for chromosome in sorted(ratio[lib][strand]):
					if memory:
						depth[lib]['1'][strand][chromosome],depth[lib]['2'][strand][chromosome]=load_varible("temp/"+lib+"/"+strand+"/depth_"+chromosome)
						ratio[lib][strand][chromosome]=load_varible("temp/"+lib+"/"+strand+"/ratio_"+chromosome)
					if gene_format:	
						if annotate == "start" and strand == '-':
							all_position=sorted(ratio[lib][strand][chromosome], reverse=True)
						elif annotate == "end" and strand == '+':
							all_position=sorted(ratio[lib][strand][chromosome], reverse=True)
						else:
							all_position=sorted(ratio[lib][strand][chromosome])
					else:
						all_position=sorted(ratio[lib][strand][chromosome])
					if gene_format:
						no_gene=False
						if len(all_gene_position_sort[strand][chromosome]) == 0:
							gene_strand=strand_check[strand_check.index(strand)-1]
							if len(all_gene_position_sort[gene_strand][chromosome]) == 0:
								no_gene=True
						else :
							gene_strand=strand
					now=0
					for position in all_position:
						enrichment=ratio[lib][strand][chromosome][position]
						if enrichment >= peak_cut_off:
							write_file.write('\t'.join([str(position),strand,chromosome])+'\t')
							if gene_format:
								if no_gene:
									write_file.write('\t'.join(['n/a','n/a','n/a','n/a'])+'\t')
								else:
									if annotate == "start":
										if gene_strand == '+':
											while True:
												start_gene=all_gene_position_sort[gene_strand][chromosome][now]
												end_gene=gene_format[gene_strand][chromosome]['annotate'][start_gene]['end']
												name_gene=gene_format[gene_strand][chromosome]['annotate'][start_gene]['name']
												feature=gene_format[gene_strand][chromosome]['annotate'][start_gene]['feature']
												if position < start_gene:
													region='Upstream'
													dist=abs(start_gene-position)
													break
												elif position <= end_gene or position == start_gene:
													region='Intra'
													dist=abs(start_gene-position)
													num_inside=1
													for start_inside in sorted(gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside']):
														end_inside=gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside'][start_inside]['end']
														if start_inside <= position <= end_inside:
															region+=' : '+','.join(map(lambda orig_string: orig_string +'-'+str(num_inside), gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside'][start_inside]['feature'].split(',')))
															break
														num_inside+=1
													break
												elif now==len(all_gene_position_sort[gene_strand][chromosome])-1:
													region='Downstream'
													dist=abs(start_gene-position)
													break
												else:
													now+=1
										else:
											now=-1
											while True:
												start_gene=all_gene_position_sort[gene_strand][chromosome][now]
												end_gene=gene_format[gene_strand][chromosome]['annotate'][start_gene]['end']
												name_gene=gene_format[gene_strand][chromosome]['annotate'][start_gene]['name']
												feature=gene_format[gene_strand][chromosome]['annotate'][start_gene]['feature']
												if position > end_gene:
													region='Upstream'
													dist=abs(end_gene-position)
													break
												elif position >= start_gene or position == end_gene:
													region='Intra'
													dist=abs(end_gene-position)
													num_inside=1
													for start_inside in sorted(gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside'], reverse=True):
														end_inside=gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside'][start_inside]['end']
														if start_inside <= position <= end_inside:
															region+=' : '+','.join(map(lambda orig_string: orig_string +'-'+str(num_inside), gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside'][start_inside]['feature'].split(',')))
															break
														num_inside+=1
													break
												elif now == -len(all_gene_position_sort[gene_strand][chromosome]):
													region='Downstream'
													dist=abs(end_gene-position)
													break
												else:
													now-=1
									else:
										if gene_strand == '-':
											while True:
												start_gene=all_gene_position_sort[gene_strand][chromosome][now]
												end_gene=gene_format[gene_strand][chromosome]['annotate'][start_gene]['end']
												name_gene=gene_format[gene_strand][chromosome]['annotate'][start_gene]['name']
												feature=gene_format[gene_strand][chromosome]['annotate'][start_gene]['feature']
												if position < start_gene:
													region='Downstream'
													dist=abs(start_gene-position)
													break
												elif position <= end_gene or position == start_gene:
													region='Intra'
													dist=abs(start_gene-position)
													num_inside=1
													for start_inside in sorted(gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside'], reverse=True):
														end_inside=gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside'][start_inside]['end']
														if start_inside <= position <= end_inside:
															region+=' : '+','.join(map(lambda orig_string: orig_string +'-'+str(num_inside), gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside'][start_inside]['feature'].split(',')))
															break
														num_inside+=1
													break
												elif now==len(all_gene_position_sort[gene_strand][chromosome])-1:
													region='Upstream'
													dist=abs(start_gene-position)
													break
												else:
													now+=1	
										else:
											now=-1
											while True:
												start_gene=all_gene_position_sort[gene_strand][chromosome][now]
												end_gene=gene_format[gene_strand][chromosome]['annotate'][start_gene]['end']
												name_gene=gene_format[gene_strand][chromosome]['annotate'][start_gene]['name']
												feature=gene_format[gene_strand][chromosome]['annotate'][start_gene]['feature']
												if position > end_gene:
													region='Downstream'
													dist=abs(end_gene-position)
													break
												elif position >= start_gene or position == end_gene:
													region='Intra'
													dist=abs(end_gene-position)
													num_inside=1
													for start_inside in sorted(gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside']):
														end_inside=gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside'][start_inside]['end']
														if start_inside <= position <= end_inside:
															region+=' : '+','.join(map(lambda orig_string: orig_string +'-'+str(num_inside), gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside'][start_inside]['feature'].split(',')))
															break
														num_inside+=1
													break
												elif now == -len(all_gene_position_sort[gene_strand][chromosome]):
													region='Upstream'
													dist=abs(end_gene-position)
													break
												else:
													now-=1
									write_file.write('\t'.join([str(dist),region,feature,name_gene,gene_strand,str(start_gene),str(end_gene)])+'\t')
							try:
								params = params_all[lib]
								arg = params[:-2]
								loc = params[-2]
								scale = params[-1] #mu
							except:
								p_value=False
							else:
								enrichment_after=boxcox_transform(enrichment,lamda)
								p_value=1-st.norm.cdf(enrichment_after,loc=loc,scale=scale, *arg)
							try:
								enrich=depth[lib]['1'][strand][chromosome][position]
							except:
								enrich=0
							try:
								unenrich=depth[lib]['2'][strand][chromosome][position]
							except:
								unenrich=0
							if distribution == 'toprank':
								write_file.write('\t'.join([str(enrich),str(unenrich),str(enrichment),'\n']))
								write_gff[strand].write('\t'.join([chromosome,lib,'.',str(position),str(position),'.',strand,'.',';'.join(['depth_1='+str(enrich),'depth_2='+str(unenrich),'ratio='+str(enrichment),'p-value='+str(p_value)])])+'\n')
							else:
								write_file.write('\t'.join([str(enrich),str(unenrich),str(enrichment),str(enrichment_after),str(p_value),'\n']))
								write_gff[strand].write('\t'.join([chromosome,lib,'.',str(position),str(position),'.',strand,'.',';'.join(['depth_1='+str(enrich),'depth_2='+str(unenrich),'ratio(before)='+str(enrichment),'ratio(after)='+str(enrichment_after),'p-value='+str(p_value)])])+'\n')
					if memory:
						depth[lib]['1'][strand][chromosome]={}
						depth[lib]['2'][strand][chromosome]={}
						ratio[lib][strand][chromosome]={}
		else:
			write_file.write(lib+" don't pass cutoff normal distribution ("+lib+"_Box-Cox_transformation.png) \n Please select -d toprank or set cutoff lower\n")
		write_file.close()

	
def write_all_position(p_value_cut_off,meta,ratio,depth,total_p_value,params_all,peak_cut_off_all,gene_format,name_out,annotate,boxcox_lamda,distribution,ratio_boxcox):
	strand_check=['+','-']
	all_gene_position_sort={}
	if gene_format:
		for strand in gene_format:
			all_gene_position_sort.update({strand:{}})
			for chromosome in gene_format[strand]:
				all_gene_position_sort[strand].update({chromosome:gene_format[strand][chromosome]['annotate'].keys()})
				all_gene_position_sort[strand][chromosome].sort()
	write_already={}
	write_file=open(meta+'_result.txt','w')
	if gene_format:
		write_file.write('position\tstrand\tchromosome\tdistance\tregion\tfeature\tf_name\tf_strand\tf_start\tf_end\t')
	else:
		write_file.write('position\tstrand\tchromosome\t')
	for lib in sorted(ratio):
		if distribution == 'normal':
			write_file.write(lib+'(depth_1)\t'+lib+'(depth_2)\t'+lib+'(ratio before boxcox)\t'+lib+'(ratio after boxcox)\t'+lib+'(p-value)\t*\t')
		else:
			write_file.write(lib+'(depth_1)\t'+lib+'(depth_2)\t'+lib+'(ratio)\t*\t')				
	if meta=="combine":
		write_file.write('combine p-value\n')
		write_gff={'+':open('_'.join(sorted(ratio))+'_combine_p-value_forward.gff','w'),'-':open('_'.join(sorted(ratio))+'_combine_p-value__reverse.gff','w')}
	else:
		write_gff={'+':open('consensus_forward.gff','w'),'-':open('consensus_reverse.gff','w')}		
		write_file.write('total datasets (pass cutoff)\n')	
	write_gff['+'].write("#gff-version 3\n")	
	write_gff['-'].write("#gff-version 3\n")
	all_lib=sorted(ratio)
	name_gff='combine_'+'_'.join(all_lib)
	for strand in total_p_value:
		write_already.update({strand:{}})
		for chromosome in sorted(total_p_value[strand]):
			write_already[strand].update({chromosome:{}})
			if memory:
				for lib2 in all_lib:
					depth[lib2]['1'][strand][chromosome],depth[lib2]['2'][strand][chromosome]=load_varible("temp/"+lib2+"/"+strand+"/depth_"+chromosome)
					ratio[lib2][strand][chromosome]=load_varible("temp/"+lib2+"/"+strand+"/ratio_"+chromosome)
			if gene_format:	
				if annotate == "start" and strand == '-':
					all_position=sorted(total_p_value[strand][chromosome], reverse=True)
				elif annotate == "end" and strand == '+':
					all_position=sorted(total_p_value[strand][chromosome], reverse=True)
				else:
					all_position=sorted(total_p_value[strand][chromosome])
			else:
				all_position=sorted(total_p_value[strand][chromosome])
			if gene_format:
				no_gene=False
				if len(all_gene_position_sort[strand][chromosome]) == 0:
					gene_strand=strand_check[strand_check.index(strand)-1]
					if len(all_gene_position_sort[gene_strand][chromosome]) == 0:
						no_gene=True
				else :
					gene_strand=strand
			now=0
			for position in all_position:
				if meta=='combine':
					combine_p=total_p_value[strand][chromosome][position]
				try:
					check=write_already[strand][chromosome][position]
				except :
					write_file.write('\t'.join([str(position),strand,chromosome])+'\t')
					if gene_format:
						if no_gene:
							write_file.write('\t'.join(['n/a','n/a','n/a','n/a'])+'\t')
						else:
							if annotate == "start":
								if gene_strand == '+':
									while True:
										start_gene=all_gene_position_sort[gene_strand][chromosome][now]
										end_gene=gene_format[gene_strand][chromosome]['annotate'][start_gene]['end']
										name_gene=gene_format[gene_strand][chromosome]['annotate'][start_gene]['name']
										feature=gene_format[gene_strand][chromosome]['annotate'][start_gene]['feature']
										if position < start_gene:
											region='Upstream'
											dist=abs(start_gene-position)
											break
										elif position <= end_gene or position == start_gene:
											region='Intra'
											dist=abs(start_gene-position)
											num_inside=1
											for start_inside in sorted(gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside']):
												end_inside=gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside'][start_inside]['end']
												if start_inside <= position <= end_inside:
													region+=' : '+','.join(map(lambda orig_string: orig_string +'-'+str(num_inside), gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside'][start_inside]['feature'].split(',')))
													break
												num_inside+=1
											break
										elif now==len(all_gene_position_sort[gene_strand][chromosome])-1:
											region='Downstream'
											dist=abs(start_gene-position)
											break
										else:
											now+=1
								else:
									now=-1
									while True:
										start_gene=all_gene_position_sort[gene_strand][chromosome][now]
										end_gene=gene_format[gene_strand][chromosome]['annotate'][start_gene]['end']
										name_gene=gene_format[gene_strand][chromosome]['annotate'][start_gene]['name']
										feature=gene_format[gene_strand][chromosome]['annotate'][start_gene]['feature']
										if position > end_gene:
											region='Upstream'
											dist=abs(end_gene-position)
											break
										elif position >= start_gene or position == end_gene:
											region='Intra'
											dist=abs(end_gene-position)
											num_inside=1
											for start_inside in sorted(gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside'], reverse=True):
												end_inside=gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside'][start_inside]['end']
												if start_inside <= position <= end_inside:
													region+=' : '+','.join(map(lambda orig_string: orig_string +'-'+str(num_inside), gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside'][start_inside]['feature'].split(',')))
													break
												num_inside+=1
											break
										elif now == -len(all_gene_position_sort[gene_strand][chromosome]):
											region='Downstream'
											dist=abs(end_gene-position)
											break
										else:
											now-=1
							else:
								if gene_strand == '-':
									while True:
										start_gene=all_gene_position_sort[gene_strand][chromosome][now]
										end_gene=gene_format[gene_strand][chromosome]['annotate'][start_gene]['end']
										name_gene=gene_format[gene_strand][chromosome]['annotate'][start_gene]['name']
										feature=gene_format[gene_strand][chromosome]['annotate'][start_gene]['feature']
										if position < start_gene:
											region='Downstream'
											dist=abs(start_gene-position)
											break
										elif position <= end_gene or position == start_gene:
											region='Intra'
											dist=abs(start_gene-position)
											num_inside=1
											for start_inside in sorted(gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside'], reverse=True):
												end_inside=gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside'][start_inside]['end']
												if start_inside <= position <= end_inside:
													region+=' : '+','.join(map(lambda orig_string: orig_string +'-'+str(num_inside), gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside'][start_inside]['feature'].split(',')))
													break
												num_inside+=1
											break
										elif now==len(all_gene_position_sort[gene_strand][chromosome])-1:
											region='Upstream'
											dist=abs(start_gene-position)
											break
										else:
											now+=1	
								else:
									now=-1
									while True:
										start_gene=all_gene_position_sort[gene_strand][chromosome][now]
										end_gene=gene_format[gene_strand][chromosome]['annotate'][start_gene]['end']
										name_gene=gene_format[gene_strand][chromosome]['annotate'][start_gene]['name']
										feature=gene_format[gene_strand][chromosome]['annotate'][start_gene]['feature']
										if position > end_gene:
											region='Downstream'
											dist=abs(end_gene-position)
											break
										elif position >= start_gene or position == end_gene:
											region='Intra'
											dist=abs(end_gene-position)
											num_inside=1
											for start_inside in sorted(gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside']):
												end_inside=gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside'][start_inside]['end']
												if start_inside <= position <= end_inside:
													region+=' : '+','.join(map(lambda orig_string: orig_string +'-'+str(num_inside), gene_format[gene_strand][chromosome]['annotate'][start_gene]['inside'][start_inside]['feature'].split(',')))
													break
												num_inside+=1
											break
										elif now == -len(all_gene_position_sort[gene_strand][chromosome]):
											region='Upstream'
											dist=abs(end_gene-position)
											break
										else:
											now-=1
							write_file.write('\t'.join([str(dist),region,feature,name_gene,gene_strand,str(start_gene),str(end_gene)])+'\t')
					pass_p=0
					for lib2 in all_lib:
						peak_cut_off=peak_cut_off_all[lib2]
						try:
							enrichment=ratio[lib2][strand][chromosome][position]
						except :
							try:
								enrich=depth[lib2]['1'][strand][chromosome][position]
								if distribution == 'normal':
									write_file.write('\t'.join([str(enrich),0,'n/a','n/a','n/a','*'])+'\t')
								else:
									write_file.write('\t'.join([str(enrich),0,'n/a','*'])+'\t')									
							except:
								if distribution == 'normal':
									try:
										unenrich=depth[lib2]['2'][strand][chromosome][position]
										write_file.write('\t'.join([0,str(unenrich),'n/a','n/a','n/a','*'])+'\t')
									except:
										write_file.write('\t'.join(['n/a','n/a','n/a','*'])+'\t')
								else:
									try:
										unenrich=depth[lib2]['2'][strand][chromosome][position]
										write_file.write('\t'.join([0,str(unenrich),'n/a','n/a','*'])+'\t')
									except:
										write_file.write('\t'.join(['n/a','n/a','n/a','*'])+'\t')
						else:
							try:
								params = params_all[lib2]
								arg = params[:-2]
								loc = params[-2]
								scale = params[-1] #mu
							except:
								p_value=False
							else:
								try:
									enrichment_after=ratio_boxcox[lib2][strand][chromosome][enrichment]
								except:
									lamda=boxcox_lamda[lib2]
									enrichment_after=boxcox_transform(enrichment,lamda)
									ratio_boxcox[lib2][strand][chromosome].update({enrichment:enrichment_after})
								p_value=1-st.norm.cdf(enrichment_after,loc=loc,scale=scale, *arg)
							try:
								enrich=depth[lib2]['1'][strand][chromosome][position]
							except:
								enrich=0
							try:
								unenrich=depth[lib2]['2'][strand][chromosome][position]
							except:
								unenrich=0
							if distribution == 'normal':
								write_file.write('\t'.join([str(enrich),str(unenrich),str(enrichment),str(enrichment_after),str(p_value)])+"\t*\t")
							else:
								write_file.write('\t'.join([str(enrich),str(unenrich),str(enrichment)])+"\t*\t")									
							if enrichment >= peak_cut_off:
								pass_p+=1
					write_already[strand][chromosome].update({position:True})
					if meta=='combine':
						write_file.write(str(combine_p)+'\n')	
						write_gff[strand].write('\t'.join([chromosome,name_gff,'.',str(position),str(position),'.',strand,'.','p-value='+str(combine_p)])+'\n')
					else:
						write_file.write(str(pass_p)+'\n')	
						write_gff[strand].write('\t'.join([chromosome,name_gff,'.',str(position),str(position),'.',strand,'.','among libraries significant='+str(pass_p)])+'\n')
			if memory:
				for lib2 in all_lib:
					depth[lib2]['1'][strand][chromosome]={}
					depth[lib2]['2'][strand][chromosome]={}
					ratio[lib2][strand][chromosome]={}					
	write_file.close()

def combine_p_value(ratio,params_all,p_value_cut_off,boxcox_lamda):
	total_p_value={'+':{},'-':{}}
	all_lib=sorted(ratio)
	ratio_boxcox={}
	for lib in ratio:
		ratio_boxcox.update({lib:{'+':{},'-':{}}})
		for strand in ratio[lib]:
			for chromosome in ratio[lib][strand]:
				ratio_boxcox[lib][strand].update({chromosome:{}})
	if memory:
		for lib in sorted(ratio):
			for strand in ratio[lib]:
				for chromosome in sorted(ratio[lib][strand]):
					ratio[lib][strand][chromosome]=load_varible("temp/"+lib+"/"+strand+"/ratio_"+chromosome)
		
	for strand in ratio[all_lib[0]]:
		for chromosome in ratio[all_lib[0]][strand]:
			total_p_value[strand].update({chromosome:{}})
			
	for lib in sorted(ratio):
		for strand in ratio[lib]:
			for chromosome in sorted(ratio[lib][strand]):
				for position in sorted(ratio[lib][strand][chromosome]):
					try:
						ckeck=total_p_value[strand][chromosome][position]
					except:
						all_p_value=[]
						for lib2 in sorted(ratio):
							try:
								enrichment=ratio[lib2][strand][chromosome][position]
							except:
								pass
							else:
								try:
									params=params_all[lib2]
									arg = params[:-2]
									loc = params[-2]
									scale = params[-1] #mu
								except:
									return {},{}
								else:
									try:
										enrichment_after=ratio_boxcox[lib2][strand][chromosome][enrichment]
									except:
										lamda=boxcox_lamda[lib2]
										enrichment_after=boxcox_transform(enrichment,lamda)
										ratio_boxcox[lib2][strand][chromosome].update({enrichment:enrichment_after})
									p_value=st.norm.cdf(enrichment_after,loc=loc,scale=scale, *arg)
								all_p_value.append(p_value)
						not_use,p_value=st.combine_pvalues(all_p_value)
						p_value=1-p_value
						if p_value <= p_value_cut_off:
							total_p_value[strand][chromosome].update({position:p_value})
	if memory:
		for lib in sorted(ratio):
			for strand in ratio[lib]:
				for chromosome in sorted(ratio[lib][strand]):
					ratio[lib][strand][chromosome]={}	
	return total_p_value,ratio_boxcox

def replicate_p_value(ratio,params_all,order,p_value_cut_off,boxcox_lamda,peak_cut_off_all):
	all_lib=sorted(ratio)
	total_p_value={'+':{},'-':{}}
	ratio_boxcox={}	
	for lib in ratio:
		ratio_boxcox.update({lib:{'+':{},'-':{}}})
		for strand in ratio[lib]:
			for chromosome in ratio[lib][strand]:
				ratio_boxcox[lib][strand].update({chromosome:{}})		
	if memory:
		for lib in sorted(ratio):
			for strand in ratio[lib]:
				for chromosome in sorted(ratio[lib][strand]):
					ratio[lib][strand][chromosome]=load_varible("temp/"+lib+"/"+strand+"/ratio_"+chromosome)
	for strand in ratio[all_lib[0]]:
		for chromosome in ratio[all_lib[0]][strand]:
			total_p_value[strand].update({chromosome:{}})
	for lib in sorted(ratio):
		for strand in ratio[lib]:
			for chromosome in sorted(ratio[lib][strand]):
				for position in sorted(ratio[lib][strand][chromosome]):
					try:
						ckeck=total_p_value[strand][chromosome][position]
					except:
						all_p_value=[]
						for lib2 in sorted(ratio):
							peak_cut_off=peak_cut_off_all[lib2]
							try:
								enrichment=ratio[lib2][strand][chromosome][position]
							except:
								pass
							else:
								try:
									params=params_all[lib2]
									arg = params[:-2]
									loc = params[-2]
									scale = params[-1] #mu
								except:
									if enrichment >= peak_cut_off:
										p_value=(1-p_value_cut_off)
									else:
										p_value=False
								else:
									try:
										enrichment_after=ratio_boxcox[lib2][strand][chromosome][enrichment]
									except:
										lamda=boxcox_lamda[lib2]
										enrichment_after=boxcox_transform(enrichment,lamda)
										ratio_boxcox[lib2][strand][chromosome].update({enrichment:enrichment_after})
									p_value=st.norm.cdf(enrichment_after,loc=loc,scale=scale, *arg)
								if p_value:
									all_p_value.append(p_value)
						all_p_value.sort()
						if order <= len(all_p_value):
							p_value=all_p_value[-order]
							p_value=1-p_value
							if p_value <= p_value_cut_off:
								total_p_value[strand][chromosome].update({position:1-all_p_value[-1]})
	if memory:
		for lib in sorted(ratio):
			for strand in ratio[lib]:
				for chromosome in sorted(ratio[lib][strand]):
					ratio[lib][strand][chromosome]={}
	return total_p_value,ratio_boxcox

def read_per_position(data1,data2,lib):
	data1=map(np.log10,data1)
	data2=map(np.log10,data2)
	density = st.gaussian_kde(data1)
	xs = np.linspace(min(data1),max(data1),100)
	density.covariance_factor = lambda : .25
	density._compute_covariance()
	plt.plot(xs,density(xs),label="Library 1",color='salmon')
	plt.fill_between(xs,0, density(xs), facecolor='salmon', alpha=0.5)
	density = st.gaussian_kde(data2)
	xs = np.linspace(min(data2),max(data2),100)
	density.covariance_factor = lambda : .25
	density._compute_covariance()
	plt.plot(xs,density(xs),label="Library 2",color='deepskyblue')
	plt.fill_between(xs,0, density(xs), facecolor='deepskyblue', alpha=0.5)
	x_axis=[[],[]]
	axes = plt.gca()
	x_axis[0]=map(int,axes.get_xticks().tolist())
	for n in x_axis[0]:
		x_axis[1].append(r'$10^'+str(n)+'$')
	plt.xticks(x_axis[0],x_axis[1])
	plt.xlabel('number of reads per position')
	plt.ylabel('Probability density')
	a=plt.title(lib)
	plt.legend(loc='upper right', shadow=True,fontsize=18)
	plt.savefig(lib+'_reads_per_Position.png', dpi=300)
	plt.clf()

# def print_out(text):
# 	global done
# 	erase=' '*(len(done)+10)
# 	sys.stdout.write("\033[F")
# 	print (erase)
# 	sys.stdout.write("\033[F")
# 	print (text+"\n\n")
# 	file_log.write(text+'\n')
# 	done=""
# def print_process(text):
# 	global done
# 	erase=' '*(len(done)+10)
# 	sys.stdout.write("\033[F")
# 	print (erase)
# 	sys.stdout.write("\033[F")
# 	print (text)
# 	done=text
def print_out(text,tell):
	global done
	if done:
		print ('Done..')
	print ("\n"+text)
	file_log.write(text+'\n')
	if tell:
		done=True
	else:
		done=False
def print_process(text,tell):
	global done
	if done:
		print ('Done..')
	print (text)
	if tell:
		done=True
	else:
		done=False
def save_varible(varible,name):
	with open(name+".pickle", 'w') as f:
		pickle.dump(varible, f)

def load_varible(name):
	with open(name+".pickle") as f:
		a = pickle.load(f)
	return a

def main():
	parser = OptionParser()
	parser.add_option("-i", "--input", dest="input", default="", type="string", 
		help="destination directory of input file")
	parser.add_option("-o", "--output", dest="output", default="ToNER", type="string",
		help="output folder name")
	parser.add_option("-r", "--reads", dest="reads", default="all", type="string",
		help="position of reads for calculate (start,all,end) ")
	parser.add_option("--none_pseudo", dest="none_pseudo", default=False, action="store_true", 
		help="none-pseudo count (default=False)" )
	parser.add_option("-t", "--total_read", dest="total_read", default=1, type="int",
		help="filter reads count with total reads (sum of both libraries)")
	parser.add_option("-e","--each_library", dest="each_library", default=False, type="int",
		help="filter reads count with minimum read depth in either library")
	parser.add_option("-p", "--p_value", dest="p_value", default=0.05, type="float",
		help="p_value cut off (default=0.05)")	
	parser.add_option("-g", "--gene",   dest="gene", default=False, type="string",
		help="add gene annotate to result [start,end] (default=False)" )
	parser.add_option("-d", "--distribution", dest="distribution", default="normal", type="string",
		help="type of distribution (normal,toprank)")
	parser.add_option("-q", "--qqplot", dest="pass_normal", default=0.9, type="float",
		help="cutoff r-squared to pass normal distribution after boxcox (default=0.9)")
	parser.add_option("-c", "--combine", dest="combine", default=False, action="store_true",
		help="meta-analysis by Fisher's combined probability test to calculate combined p-values")
	parser.add_option("-s", "--consensus", dest="consensus", default=False, type="int",
		help="called significant among at least a minimum number of replicates specified ")
	parser.add_option("--less_memory", dest="memory", default=False, action="store_true",
		help="In case of memory problem, this option will use less memory than default which that are time consuming")
	parser.add_option("--history_data", dest="history_data", default=False, action="store_true",
		help="Added all output file from each process of ToNER program")
	
	(options, args) = parser.parse_args()
	
	
	global file_log,memory,done
	depth={}
	ratio={}
	params_all={}
	boxcox_lamda={}
	all_lib=[]
	all_file=[]
	percentile_lib={}
	peak_cut_off_all={}
	chromosome_lib={}
	if not options.input:
		parser.print_help()
		sys.exit("Please put input directory in -i option")
	input_folder=options.input
	position_depth=options.reads
	p_value_cut_off=options.p_value
	annotate=options.gene
	combine=options.combine
	consensus=options.consensus
	output=options.output
	pass_normal=options.pass_normal
	distribution=options.distribution
	memory=options.memory
	save_data=options.history_data
	total_read=options.total_read
	each_library=options.each_library
	if options.none_pseudo:
		pseudo=False
	else:
		pseudo=True
	if each_library:
		filter_type="each"
		reads_count=each_library
	else:
		filter_type="total"
		reads_count=total_read
	if combine:
		meta="combine"
	elif consensus:
		meta="consensus"
		replicate=consensus
	else:
		meta=False
	if distribution not in ['normal','toprank']:
		sys.exit('!error : -d,--distribution  must be "normal" or "toprank"')
	if annotate and annotate not in ['start', 'end']:
		sys.exit('!error : -g,--gene must be "start" or "end"')
	if position_depth not in ['start','all','end']:
		sys.exit('!error : -r,--reads must be "start", "all" or "end"')
	if input_folder[-1] != '/':
		input_folder=input_folder+'/'
	if os.path.isdir(input_folder):
		all_file+=glob.glob(os.path.join(input_folder, '*.bam'))
	else:
		sys.exit('!error : No such directory ('+input_folder+')')
	all_file.sort()
	for x in range(len(all_file)):
		all_file[x]=all_file[x].split("/")[-1]
	if len(all_file) == 0:
		sys.exit("!error : No such bam file in this folder("+input_folder+")")
	else:
		file1=[]
		file2=[]
		for file in all_file:
			if '_1.bam' in file or '_2.bam' in file:
				if '_1.bam' in file:
					file1.append(file.replace('_1.bam',''))
				elif '_2.bam' in file:
					file2.append(file.replace('_2.bam',''))
			else :
				sys.exit('!error : Invalid file name ('+file+'), File name must end with "_1.bam or _2.bam" (exp. ecoli_1.bam or ecoli_2.bam)')
		if len(set(file1)-set(file2)) != 0:
			for file in set(file1)-set(file2):
				print("No such file "+file+"_2.bam")
			sys.exit()
		elif len(set(file2)-set(file1)) != 0:
			for file in set(file2)-set(file1):
				print("No such file "+file+"_1.bam")
			sys.exit()
	if consensus and len(all_file)/2 < consensus:
		sys.exit('!error : number of datasets pass p-value more than input datasets')
	elif consensus and consensus < 1:
		sys.exit('!error : number of datasets pass p-value must more than 0')
	if pass_normal > 1:
		pass_normal=1
	# ********************start here****************
	time_start=str(datetime.now())[:-7].replace(' ','_').replace(':','_')
	
	if output != 'ToNER':
		make_directory(input_folder+output,save_data)
	else:
		make_directory(input_folder+output+'_'+time_start,save_data)
	
	time_start=time.time()	
	
	file_log=open('log.txt','w')
	done=''
	print_out("ToNER 1.0\n",False)
	print_out("\ncommand-line : "+' '.join(map(str,sys.argv))+'\n',False)	
	print_out("\n--------------------------",False)
	print_out("|"+(str(datetime.now())[:-7])+"|",False)
	print_out("--------------------------\n",False)
	print_out("input_folder :\t"+input_folder,False)
	print_out("output folder :\t" +os.getcwd(),False)
	print_out(" " +str(int(len(all_file)/2))+ " datasets :",False)
	for file in file1:
		print_out('\t- '+file,False)
	if meta:
		if combine:
			print_out("meta-analysis :\tcombine p-value",False)
		elif consensus:
			print_out("meta-analysis :\tconsensus "+str(consensus)+ " datasets",False)
	print_out("p-value cut off :\t"+ str(p_value_cut_off),False)
	print_out("position of reads :\t"+position_depth,False)
	print_out("type of filter :\t"+str(filter_type),False)
	print_out("minimum of reads count :\t"+str(reads_count),False)
	if pseudo: 
		print_out("pseudo count :\tYes",False)
	else: 
		print_out("pseudo count :\tNo",False)
	if distribution == "normal":
		print_out("distribution :\tnormal distribution",False)
		print_out("r-squared cutoff for pass normal distribution after Box-Cox transformation :\t"+str(pass_normal),False)
	else:
		print_out("distribution :\ttoprank",False)
	
	if annotate : 
		print_out("put annotate gene :\tgene "+annotate,False)
	else: 
		print_out("put annotate gene :\tNo",False)
	if memory:
		print_out("less memory :\tYes (use long time)",False)
	if save_data:
		print_out('save raw data in folder "history_data"',False)
	print_out("\n****************************\n\n",False)
	
	if annotate:
		print_process("import gene annotate...",True)
		gene_format,num_gene=input_gene_gff(input_folder)
		if not num_gene:
			annotate=False

	for file1 in all_file:
		if '_1.' in file1:
			lib=file1[:-6]
			if lib not in all_lib:
				print_out("\n"+lib+" :\n",False) 
				all_lib.append(lib)
				depth[lib]={}
				ratio[lib]={}
				file2=file1.replace('_1.','_2.')
				if memory:
					if not os.path.isdir("temp/"+lib):
						os.makedirs("temp/"+lib)
					if not os.path.isdir("temp/"+lib+"/+"):
						os.makedirs("temp/"+lib+"/+")
						os.makedirs("temp/"+lib+"/-")
				depth[lib],chromosome_lib[lib] = data_process(file1,file2,lib,position_depth)
				print_process("calculate ratio "+lib+"...",True)
				if pseudo:
					both,only_1,only_2=check_pseudo(depth[lib]['1'],depth[lib]['2'],reads_count,filter_type)
					depth[lib]['total_reads_1']+=only_2
					depth[lib]['total_reads_2']+=only_1
								
				ratio[lib],both,only_1,only_2=find_ratio(
									depth[lib]['1'],depth[lib]['total_reads_1'],
									depth[lib]['2'],depth[lib]['total_reads_2'],
									pseudo,
									lib,
									save_data,
									reads_count,
									filter_type
									)
				pie_plot([both,only_1,only_2],['both','only_1','only_2'],lib+'_type_of_reads',lib)
				# read_per_position(data1,data2,lib)
				if distribution == "normal":
					print_process("transform data with boxcox method "+ lib +"...",False)
					boxcox_lamda[lib],params_all[lib],peak_cut_off_all[lib],peak_after,r2,normal = boxcox_data_transform(ratio[lib],lib,p_value_cut_off,pass_normal)
					print_out("\t- r-squared QQ-plot = "+str(r2),False)
					if normal: 
						mu, std = params_all[lib]
						print_out("\t- parameters :",False)
						print_out("\t\tmu = "+str(mu),False)
						print_out("\t\tstd = "+str(std),False)
						print_out("\t\tlambda = "+str(boxcox_lamda[lib]),False)
						print_out("\t\tpeak cutoff = "+str(peak_cut_off_all[lib])+" ("+str(peak_after)+")",False)
						save_history_data(ratio[lib],boxcox_lamda,lib,save_data,True)
					else : 
						del params_all[lib],boxcox_lamda[lib]
						if save_data:
							save_history_data(ratio[lib],boxcox_lamda,lib,save_data,False)
						print_out("\tdon't pass cutoff normal distribution ("+lib+"_Box-Cox_transformation.png) Please select -d toprank or set cutoff lower",False)
				else:
					if save_data:
						save_history_data(ratio[lib],boxcox_lamda,lib,save_data,False)
					peak_cut_off_all[lib]=find_empirical(ratio[lib],(1-p_value_cut_off)*100,lib)
					print_out("\tpeak cutoff = "+str(peak_cut_off_all[lib]),False)
	if annotate:
		for lib in chromosome_lib.keys():
			for chromosome in chromosome_lib[lib].keys():
				length=chromosome_lib[lib][chromosome]
				if chromosome in gene_format['+']:
					num=gene_format['+'][chromosome]['length']
					if num != length:
						annotate=False
						print_out("\ngene annotate :\n\tcouldn't put gene annotate to output file because in "+lib+" length chromosome '"+chromosome+"' not match with gene annotate in gff file (gff="+str(num)+", "+lib+"="+str(length)+")",False)
						break
				else:
					annotate=False
					print_out("\ngene annotate :\n\tcouldn't put gene annotate to output file because couldn't find chromosome '"+chromosome+"' in gene annotate",False)
					print_out("\tall chromosome of gff:",False)
					for chromosome_gff in sorted(gene_format['+']):
						print_out("\t\t"+chromosome_gff,False)
					print_out("\tall chromosome of : "+lib+" alignment",False)
					for chromosome_gff in sorted(chromosome_lib[lib].keys()):
						print_out("\t\t"+chromosome_gff,False)
					break
			else:
				continue
			break
	if annotate:
		write_defult_result(ratio,depth,p_value_cut_off,params_all,peak_cut_off_all,gene_format,annotate,boxcox_lamda,distribution)
	else:
		write_defult_result(ratio,depth,p_value_cut_off,params_all,peak_cut_off_all,False,annotate,boxcox_lamda,distribution)
	total_p_value={}
	if meta=='combine' and distribution == 'toprank':
		print_out("\nmeta-analysis :\n\tcouldn't combine p-value for toprank",False)
	elif meta and len(all_file)> 2:
		if meta == 'combine':
			can=True
			for lib in all_lib:
				if lib not in boxcox_lamda:
					can=False
					break
			if can:
				print_process("combine p-value...",True)
				total_p_value,ratio_boxcox=combine_p_value(ratio,params_all,p_value_cut_off,boxcox_lamda)
			else:
				print_out("\nmeta-analysis :\n\tcouldn't combine p-value because some of dataset don't pass cutoff normal distribution",False)
		elif meta == 'consensus':
			can=True
			if distribution == 'normal':
				for lib in all_lib:
					if lib not in boxcox_lamda:
						can=False
						break
			if can:
				print_process("consensus datasets..",True)
				total_p_value,ratio_boxcox=replicate_p_value(ratio,params_all,replicate,p_value_cut_off,boxcox_lamda,peak_cut_off_all)
			else:
				print_out("\nmeta-analysis :\n\tcouldn't consensus because some of dataset don't pass cutoff normal distribution",False)
		if total_p_value:
			if annotate:
				print_process("write result "+meta+"...",True)
				write_all_position(p_value_cut_off,meta,ratio,depth,total_p_value,params_all,peak_cut_off_all,gene_format,output,annotate,boxcox_lamda,distribution,ratio_boxcox)
			else:
				print_process("write result "+meta+"...",True)
				write_all_position(p_value_cut_off,meta,ratio,depth,total_p_value,params_all,peak_cut_off_all,False,output,annotate,boxcox_lamda,distribution,ratio_boxcox)
	elif meta and len(all_file) ==2 :
		print_out("\nmeta-analysis :\n\tprogram couldn't meta analysis with 1 datasets",False)
	time_use=(time.time()-time_start)
	print_out("\nTime " + str(time_use) + " secs.",False)
	print_out("\n--------------------------",False)
	print_out("|"+(str(datetime.now())[:-7])+"|",False)
	print_out("--------------------------",False)
	file_log.close()
	shutil.rmtree("temp")
	
if __name__ == '__main__' :
	main()
	

