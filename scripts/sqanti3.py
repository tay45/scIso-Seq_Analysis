#!/usr/bin/python
import os
import os.path
import subprocess
import fileinput
import argparse
import sys
import time
import gzip
import shutil



#Check the input files (no parser) 
def input_file2 (path):
	if os.path.isfile(path):
		return
	else:
		print("Doesn't exist. Please check the file location")
		sys.exit()



#Parameters
input_gtf = "dedup.5merge.collapsed.gtf"
abundance = "dedup.5merge.collapsed.abundance.txt"
out_name = "dedup.5merge.collapsed"

#Add annotation file
annot_gtf = ""
annot_gtf = input("Enter the path of the 'annotation_gtf'.: ")
input_file2(annot_gtf)

#Add reference
ref_fasta = ""
ref_fasta = input("Enter the path of the 'reference_fasta'.: ")
input_file2(ref_fasta)

#Add cage_peak file
cage_peak = ""
cage_peak = input("Enter the path of the 'CAGE_Peak_bed'.: ")
input_file2(cage_peak)

#Add polya_peak file
polya_peak = ""
polya_peak = input("Enter the path of the 'polyA_peak_bed'.: ")
input_file2(polya_peak)

#Add polya_motif file
polya_motif = ""
polya_motif = input("Enter the path of the 'polyA_motif_txt'.: ")
input_file2(polya_motif)

#Add tappAS-annotation file
tappAS = ""
tappAS = input("Enter the path of the 'tappAS-annotation_gff3'.: ")
input_file2(tappAS)

#Add short_coverage file
coverage = ""
coverage = input("Enter the path of the 'intropolis junction_bed'.: ")
input_file2(coverage)

#CPUs
cpus = ""
cpus = input("Enter the number of the cpus to use: ")



#Convert gff to gtf
input_file2("dedup.5merge.collapsed.gff")
current_path = os.path.abspath(os.getcwd())
print(current_path)
collapsed_gff = os.path.join(current_path, "dedup.5merge.collapsed.gff")
out = os.path.join(current_path, "dedup.5merge.collapsed.gtf")
convert_gff = "gffread -T" + " " + collapsed_gff + " -o" + " " + out
os.system(convert_gff)
subprocess.call("ls dedup.*", shell = True)
time.sleep(2)



#Change the input file name
#subprocess.call("mv collapsed.abundance.txt collapsed.abundance.tsv", shell = True)
#input_file2("collapsed.abundance.tsv")



#Short Reads
while True:
	answer = input("Do you want to provide Short-Read fastq files (fofn)? (yes/no):")
	if answer.lower().startswith("n"):
		print("***Running SQANTI***")
		time.sleep(2)
		break
	elif answer.lower().startswith("y"):
		short_reads = ""
		short_reads = input("Enter the path of the 'short_reads.fastq (fofn)'.: ")
		input_file2(short_reads)
		sqanti3_qc_shortRead =  "sqanti3_qc.py " + input_gtf + " " + annot_gtf + " " + ref_fasta + " --cage_peak " + cage_peak + " --polyA_peak " + polya_peak + " --polyA_motif_list " + polya_motif + " -o " + out_name + " -d " + "sqanti3" + " --fl_count " + abundance + " --short_reads " + short_reads + " --cpus " + str(cpus) + " --genename" + " --isoAnnotLite" + " --gff3 " + tappAS + " --report both"
		os.system(sqanti3_qc_shortRead)
		#subprocess.call("mv collapsed.abundance.tsv collapsed.abundance.txt", shell = True)
		
		
		
		#Run gtfToGenePred
		print("***Run gtfToGenePred***")
		time.sleep(2)
		input_file2("dedup.5merge.collapsed.gtf")
		gtfToGenePred = "gtfToGenePred" + " " + "dedup.5merge.collapsed.gtf" + " " + "dedup.5merge.collapsed.genePred"
		os.system(gtfToGenePred)



		#Run genePredToBed
		print("***Run genePredToBed***")
		time.sleep(2)
		input_file2("dedup.5merge.collapsed.genePred")
		genePredToBed = "genePredToBed" + " " + "dedup.5merge.collapsed.genePred" + " " + "dedup.5merge.collapsed.bed12"
		os.system(genePredToBed)



		#Run color_bed12
		print("***Run color_bed12_post_sqanti.py***")
		time.sleep(2)
		#Current working directory
		current_path = os.path.abspath(os.getcwd())
		print(current_path)
		input_file2(current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_classification.txt")
		color_bed12 = "color_bed12_post_sqanti.py" + " " + current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_classification.txt" + " " + "dedup.5merge.collapsed.bed12" + " " + "dedup.5merge.collapsed.colored"
		os.system(color_bed12)
		subprocess.call("ls dedup.*", shell = True)
		time.sleep(2)
		
		
		
		#Run IsoAnnot
		print("***Run IsoAnnot***")
		time.sleep(2)
		#Current working directory
		current_path = os.path.abspath(os.getcwd())
		print(current_path)
		#Check file locations
		input_file2(current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_corrected.gtf")
		input_file2(current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_classification.txt")
		input_file2(current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_junctions.txt")
		correct = current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_corrected.gtf"
		classification = current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_classification.txt"
		junctions = current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_junctions.txt"
		#Run IsoAnnot.py
		isoannot = "python3 IsoAnnotLite_v2.7.0_SQ3.py" + " " + correct + " " + classification + " " + junctions + " " + "-gff3" + " " + tappAS + " " + "-novel -o" + " " + "isoannot" + " " + "-stdout" + " " + "isoannot.summaryResults"
		os.system(isoannot)
		subprocess.call("ls isoannot.*", shell = True)
		time.sleep(2)
		


#Run SQANTI3
sqanti3_qc = "sqanti3_qc.py " + input_gtf + " " + annot_gtf + " " + ref_fasta + " --cage_peak " + cage_peak + " --polyA_peak " + polya_peak + " --polyA_motif_list " + polya_motif + " -o " + out_name + " -d " + "sqanti3" + " --fl_count " + abundance + " --cpus " + str(cpus) + " --genename" + " --isoAnnotLite" + " --gff3 " + tappAS + " -c " + coverage + " --report both"
subprocess.run(sqanti3_qc, shell=True)
#subprocess.call("mv collapsed.abundance.tsv collapsed.abundance.txt", shell = True)



#Run gtfToGenePred
print("***Run gtfToGenePred***")
time.sleep(2)
input_file2("dedup.5merge.collapsed.gtf")
gtfToGenePred = "gtfToGenePred" + " " + "dedup.5merge.collapsed.gtf" + " " + "dedup.5merge.collapsed.genePred"
os.system(gtfToGenePred)



#Run genePredToBed
print("***Run genePredToBed***")
time.sleep(2)
input_file2("dedup.5merge.collapsed.genePred")
genePredToBed = "genePredToBed" + " " + "dedup.5merge.collapsed.genePred" + " " + "dedup.5merge.collapsed.bed12"
os.system(genePredToBed)



#Run color_bed12
print("***Run color_bed12_post_sqanti.py***")
time.sleep(2)
#Current working directory
current_path = os.path.abspath(os.getcwd())
print(current_path)
input_file2(current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_classification.txt")
color_bed12 = "color_bed12_post_sqanti.py" + " " + current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_classification.txt" + " " + "dedup.5merge.collapsed.bed12" + " " + "dedup.5merge.collapsed.colored"
os.system(color_bed12)
subprocess.call("ls dedup.*", shell = True)
time.sleep(2)
		
		
		
#Run IsoAnnot
print("***Run IsoAnnot***")
time.sleep(2)
#Current working directory
current_path = os.path.abspath(os.getcwd())
print(current_path)
#Check file locations
input_file2(current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_corrected.gtf")
input_file2(current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_classification.txt")
input_file2(current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_junctions.txt")
correct_gtf = current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_corrected.gtf"
classification = current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_classification.txt"
junctions = current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_junctions.txt"
#Run IsoAnnot.py
isoannot = "python3 IsoAnnotLite_v2.7.0_SQ3.py" + " " + correct_gtf + " " + classification + " " + junctions + " " + "-gff3" + " " + tappAS + " " + "-novel -o" + " " + "isoannot" + " " + "-stdout" + " " + "isoannot.summaryResults"
os.system(isoannot)
subprocess.call("ls isoannot.*", shell = True)
time.sleep(2)



#Filter Artifacts
print("***Filter Artifacts***")
time.sleep(2)
#Current working directory
current_path = os.path.abspath(os.getcwd())
print(current_path)
#Check file locations
input_file2(current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_classification.txt")
input_file2(current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_corrected.fasta")
input_file2("dedup.5merge.collapsed.gff")
classification = current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_classification.txt"
correct_fasta = current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_corrected.fasta"
gff = current_path + "/" + "dedup.5merge.collapsed.gff"
#Run IsoAnnot.py
rulefilter = "sqanti3_RulesFilter.py" + " " + classification + " " + correct_fasta + " " + gff
os.system(rulefilter)
subprocess.call("ls" + " " + current_path + "/sqanti3/" + "dedup.5merge.collapsed_classification.filtered_lite.*", shell = True)
time.sleep(2)
  