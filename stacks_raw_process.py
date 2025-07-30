#!/usr/bin/env python

#This script can be used to clean and de-multiplex RADseq data using Stacks and cutadapt (both required). It applies the clone_filter function to remove PCR duplicates. Then it uses cutadapt to trim adapters (usually not necessary) and to trim the reads to certain lengths and remove short reads. Finally, it uses process_radtags to de-multiplex the RADseq libraries. By switching options on and off only specific steps can be re-run, but generally the steps build on each other and everything should be run.

#Input: Path to directory containing the fastq files

import sys
import os 
from glob import glob


clone_filt=False #switch on or off as needed
trim=False
process_radtags=True


path=sys.argv[1]
os.chdir(path)

filelist=glob("*.fq.gz") #adjust as needed

filelist.sort() #sort the list by name so that the second read file always follows directly the first read file (maybe check if this is true for your files) 

try:
	os.mkdir("duprem/")
except:
	a=None	
try:
	os.mkdir("samples/")
except:
	a=None
try:
	os.mkdir("trimmed/")
except:
	a=None		

c=0
cmd_list=[]
if clone_filt:
	for i in range(0,len(filelist),2): #the loop assumes that the fastq files are ordered so that the second read file always follows directly the first read file (should be the normal case)
		c+=1
		file1=filelist[i]
		file2=filelist[i+1]
		clone_cmd=" clone_filter -1 {0} -2 {1} -o duprem/ -i gzfastq --oligo_len_1 5 --inline_null 2> duprem/clonfilt_log{2}".format(file1,file2,c)
		print clone_cmd
		os.system(clone_cmd)

filelist=glob("duprem/*.fq.gz")
filelist.sort()	
c=0
if trim:
	for i in range(0,len(filelist),2):
		c+=1
		file1=filelist[i]
		name1=file1.split("/")[-1]
		file2=filelist[i+1]
		name2=file2.split("/")[-1]
		trim_cmd="cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -q 30 -m 80 --trim-n -o trimmed/{0} -p trimmed/{1} {2} {3}".format(name1,name2,file1,file2)
		print trim_cmd
		os.system(trim_cmd)
		
filelist=glob("trimmed/*.fq.gz")
filelist.sort()	
c=0
if process_radtags:
	for i in range(0,len(filelist),2):
		c+=1
		file1=filelist[i]
		file2=filelist[i+1]
		process_cmd="process_radtags -1 {0} -2 {1} -e kpnI --inline_null --barcode_dist_1 1 -c -q -r -b barcodes_pool{2} -o samples/ -i gzfastq -y fastq".format(file1,file2,c) #here you may have to adjust the path/name of the barcode file; it is assumed that the pool numbers are in the same order as your fastq file (which is usually the case)
		print process_cmd
		os.system(process_cmd)
		os.system("mv samples/process_radtags.trimmed.log samples/process_radtags{0}.log".format(c))



