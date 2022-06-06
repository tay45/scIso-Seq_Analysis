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
		#sys.exit()



def read_stat():
	#Generate '.read_stat.txt'
	print("***Selecting gene loci with sufficient coverage for phasing...***")
	time.sleep(2)
	
	#Current working directory
	current_path = os.path.abspath(os.getcwd())
	print(current_path)
	
	#Check file locations
	input_file2("dedup.5merge.collapsed.group.txt")
	
	#Variables of input files
	group = os.path.join(current_path, "dedup.5merge.collapsed.group.txt")
	out = os.path.join(current_path, "dedup.5merge.collapsed.read_stat.txt")
	
	#Run 'convert_group_to_read_stat_file.py'
	stat = "python3 /opt/cDNA_Cupcake/cDNA_Cupcake/phasing/utils/convert_group_to_read_stat_file.py" + " " + group + " " + out
	os.system(stat)
	subprocess.call("ls dedup.*", shell = True)
	time.sleep(2)



def Selecting_gene_loci():
	#Selecting gene loci with sufficient coverage for phasing
	#Enter reference genome information
	ref = ""
	ref = input("Enter the path of the 'genome reference (.fa)': ")
	input_file2(ref)
	
	#Current working directory
	current_path = os.path.abspath(os.getcwd())
	print(current_path)
	
	#Check file locations
	input_file2("dedup.fasta")
	input_file2(current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_classification.filtered_lite.gtf")
	input_file2("dedup.5merge.collapsed.read_stat.txt")
	
	#Variables of input files
	fasta = os.path.join(current_path, "dedup.fasta")
	lite = current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_classification.filtered_lite.gtf"
	stat = os.path.join(current_path, "dedup.5merge.collapsed.read_stat.txt")
	
	#Run 'select_loci_to_phase.py'
	loci = "python3 /opt/cDNA_Cupcake/cDNA_Cupcake/phasing/utils/select_loci_to_phase.py" + " " + ref + " " + fasta + " " + lite + " " + stat + " -c 40"
	os.system(loci)
	subprocess.call("ls" + " " + current_path + "/" + "by_loci", shell = True)
	print("Outputs are in ../by_loci")
	time.sleep(2)
	


def Generate_and_run_IsoPhase():
	#Generate and run IsoPhase commands
	print("***Generate and run IsoPhase commands...***")
	time.sleep(2)
	
	current_path = os.path.abspath(os.getcwd())
	current_path = current_path + "/by_loci/"
	print(current_path)
	time.sleep(2)
	
	phase = ""				
	phase = input("Enter the name of a by_loci/PB.X_sizeN subdirectory:")
	print(phase)
	time.sleep(2)
	
	current_path = current_path + phase + "/"
	print(current_path)
	time.sleep(2)
	
	#Copy run_phasing_in_dir.sh
	subprocess.call("cp /opt/cDNA_Cupcake/cDNA_Cupcake/phasing/utils/run_phasing_in_dir.sh" + " " + current_path, shell = True)
	
	input_file2(current_path + "/run_phasing_in_dir.sh")
	
	#Generate a per-locus bash script
	subprocess.call("ls -1d " + current_path + "|xargs -n1 -i echo" + " " + '"bash run_phasing_in_dir.sh {}"' + " " + ">" + current_path + "/cmd", shell = True)
	subprocess.call("chmod u+x " + current_path + "/cmd", shell =  True)
	os.chdir(current_path)
	subprocess.call("bash cmd", shell = True)
	
	#Run generated a per-locus bash script
	subprocess.call("ls -1d " + current_path + "|xargs -n1 -i echo" + " " + '"cd {}; bash run.sh; cd ../../"' + " " + ">" + current_path + "cmd2", shell = True)
	subprocess.call("chmod u+x " + current_path + "cmd2", shell =  True)
	subprocess.call("bash cmd2", shell = True)
	
	subprocess.call("ls " + current_path + "*.vcf", shell = True)



def Tagging_BAM():
	#Tagging BAM files with phasing information
	print("***Tagging BAM files with phasing information...***")
	time.sleep(2)
	
	#Align to Genome
	print("***Align to Genome...***")
	time.sleep(2)
	
	current_path = os.path.abspath(os.getcwd())
	current_path = current_path + "/"
	print(current_path)
	time.sleep(2)
		
	#Check the path of input file
	input_file2(current_path + "ccs.fasta")
	
	#Add reference file
	reference = ""
	reference = input("Enter the path of the 'genome reference (.fa)': ")
	input_file2(reference)
	fasta = os.path.join(current_path, "ccs.fasta")
	out = os.path.join(current_path, "ccs.hg38.sorted.bam")
	
	#Align the input (dedup reads) to genome to get a BAM file; use GMAP or minimap2
	align = "minimap2 -ax splice -uf --secondary=no -C5" + " " + reference + " " + fasta + "|samtools view -bS|samtools sort > " + out
	os.system(align)
	subprocess.call("ls " + current_path +  "ccs.hg38.sorted.bam", shell = True)
	time.sleep(2)
	
	#Tag the BAM file
	print("***Tag the BAM file...***")
	time.sleep(2)
	
	#Check the path of input file
	input_file2(current_path + "ccs.hg38.sorted.bam")
	input_file2(current_path + "phased.partial.cleaned.hap_info.txt")
	bam = os.path.join(current_path, "ccs.hg38.sorted.bam")
	hap = os.path.join(current_path, "phased.partial.cleaned.hap_info.txt")
	out = os.path.join(current_path, "ccs.hg38.sorted.tagged.bam")
	
	#Run tag_bam_post_phasing.py
	tag = "python3 /opt/cDNA_Cupcake/cDNA_Cupcake/phasing/utils/tag_bam_post_phasing.py" + " " + bam + " " + hap + " " + out
	os.system(tag)
	
	#Run samtools
	samtools = "samtools index " + current_path + "ccs.hg38.sorted.tagged.bam"
	os.system(samtools)
	subprocess.call("ls " + current_path + "ccs.hg38.sorted.*", shell = True)
	time.sleep(2)



def Tagging_BAM_2():
	#Tagging BAM files with phasing information
	print("***Tagging BAM files with phasing information...***")
	time.sleep(2)
	
	#Align to Genome
	print("***Align to Genome...***")
	time.sleep(2)
	
	current_path = os.path.abspath(os.getcwd())
	current_path = current_path + "/"
	print(current_path)
	time.sleep(2)
	
	#Check the path of input file
	input_file2(current_path + "ccs.fasta")
	
	#Add reference file
	reference = ""
	reference = input("Enter the path of the 'genome reference (.fa)': ")
	input_file2(reference)
	fasta = os.path.join(current_path, "ccs.fasta")
	out = os.path.join(current_path, "ccs.hg38.sorted.bam")
	
	#Align the input (dedup reads) to genome to get a BAM file; use GMAP or minimap2
	align = "minimap2 -ax splice -uf --secondary=no -C5" + " " + reference + " " + fasta + "|samtools view -bS|samtools sort > " + out
	os.system(align)
	subprocess.call("ls " + current_path +  "ccs.hg38.sorted.bam", shell = True)
	time.sleep(2)
	
	#Tag the BAM file
	print("***Tag the BAM file...***")
	time.sleep(2)
	#Check the path of input file
	input_file2(current_path + "ccs.hg38.sorted.bam")
	input_file2(current_path + "phased.partial.cleaned.hap_info.txt")
	bam = os.path.join(current_path, "ccs.hg38.sorted.bam")
	hap = os.path.join(current_path, "phased.partial.cleaned.hap_info.txt")
	out = os.path.join(current_path, "ccs.hg38.sorted.tagged.bam")
	
	#Run tag_bam_post_phasing.py
	tag = "python3 /opt/cDNA_Cupcake/cDNA_Cupcake/phasing/utils/tag_bam_post_phasing.py" + " " + bam + " " + hap + " " + out + " --celltype " +	cell_type
	os.system(tag)
	
	#Run samtools
	samtools = "samtools index " + current_path + "ccs.hg38.sorted.tagged.bam"
	os.system(samtools)
	subprocess.call("ls " + current_path + "ccs.hg38.sorted.*", shell = True)
	time.sleep(2)



def summary():
	#Summarizing IsoPhase output
	#subprocess.call("cd ..", shell = True)
	#subprocess.call("cd ..", shell = True)

	#Current working directory
	current_path = os.path.abspath(os.getcwd())
	print(current_path)
	time.sleep(2)

	print("***Summarizing IsoPhase output...***")
	time.sleep(2)
	summary = "python3 /opt/cDNA_Cupcake/cDNA_Cupcake/phasing/utils/summarize_byloci_results.py"
	os.system(summary)
	subprocess.call("ls summarized.isophase_results.txt", shell = True)
	time.sleep(2)

	#Collect all the VCFs from each separate directory into a single file
	print("***collect all the VCFs from each separate directory into a single file...***")
	time.sleep(2)
	collect = "python3 /opt/cDNA_Cupcake/cDNA_Cupcake/phasing/utils/collect_all_vcf.py"
	os.system(collect)
	subprocess.call("ls IsoSeq_IsoPhase.vcf", shell = True)
	time.sleep(2)


#Advanced: linking molecule info to cell type, running IsoPhase, tagging BAMs
while True:
	answer = input("Do you have cell barcode-to-cell type information? (yes/no):")
	if answer.lower().startswith("n"):
		break
	elif answer.lower().startswith("y"):
		#Linking molecule info to cell type
		print("***Linking molecule info to cell types...***")
		time.sleep(2)
		
		#Current working directory
		current_path = os.path.abspath(os.getcwd())
		print(current_path)
		
		#Enter cell type information
		cell_type = ""
		cell_tyle = input("Enter the path of the 'bc_info.csv': ")
		input_file2(cell_type)
		
		#Check file locations
		input_file2("dedup.annotated.csv")
		
		#Variables of input files
		annot = os.path.join(current_path, "dedup.annotated.csv")
		out = os.path.join(current_path, "dedup.annotated.celltype.csv")
		
		#Run 'link_molecule_to_celltype.py' 
		link = "python3 /opt/cDNA_Cupcake/cDNA_Cupcake/singlecell/link_molecule_to_celltype.py" + " " + cell_tyle + " " + annot + " " + out
		os.system(link)
		subprocess.call("ls dedup.*", shell = True)
		time.sleep(2)
		
		read_stat()
		
		Selecting_gene_loci()
		
		while True:
			answer = input("Do you have more loci to analyze? (yes/no):")
			if answer.lower().startswith("n"):
				summary()
				print("***The running has been completed...***")
				break
			elif answer.lower().startswith("y"):
				Generate_and_run_IsoPhase()
				Tagging_BAM_2()



read_stat()

Selecting_gene_loci()

while True:
	answer = input("Do you have more loci to analyze? (yes/no):")
	if answer.lower().startswith("n"):
		summary()
		print("***The running has been completed...***")
		break
	elif answer.lower().startswith("y"):
		Generate_and_run_IsoPhase()
		Tagging_BAM()
