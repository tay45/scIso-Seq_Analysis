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


#Load modules
#module load smrtlink/11.0
#module load cDNA_Cupcake
#module load SQANTI3
#module load gffread/0.12.7



#Check the path of input files
def input_path (path):
	if os.path.isfile(path):
		return path
	else:
		parser.error("Doesn't exist the %s!" % path)
	
		
#Check the input file
def input_file (parser, arg):
	if os.path.isfile(arg):
		return open(arg, 'r')
	else:
		parser.error("Doesn't exist the %s!" % arg)


#Check the input files (no parser) 
def input_file2 (path):
	if os.path.isfile(path):
		return
	else:
		print("Doesn't exist. Please check the file location")
		sys.exit()



#Commend line arguments in lima
parser = argparse.ArgumentParser(description = "Running Iso-Seq pipeline")
parser.add_argument("-c", dest="ccs", required=True, help='path of input ccs.bam', metavar="CCS File", type=input_path)
parser.add_argument("-p", dest="primers", required=True, help='path of input primers.fasta', metavar="Primer File", type=input_path)
parser.add_argument("-o", dest="output", required=True, help='output file (fl.bam)', metavar="fl.bam")
parser.add_argument("-pf", dest="primers_file", required=False, help='input primers.fasta', metavar="Primer File", type=lambda x: input_file(parser, x))
args = parser.parse_args()



#Confirm the setting of SMRT_Root
while True:
	answer = input("Have you loaded the required modules before running the pipeline? (yes/no):")
	if answer.lower().startswith("y"):
		print("***Start the pipeline.***")
		time.sleep(2)
		break
	elif answer.lower().startswith("n"):
		print("***Please load 'smrtlink/10.1', 'cDNA_Cupcake', 'SQANTI3', and 'gffread/0.12.7' using 'module load' in Apollo.***")
		time.sleep(2)
		sys.exit()



#Deconcatenate reads for the HIT-scISOseq
while True:
	answer = input("Do you need to deconcatenate reads into their original cDNA sequences? (yes/no):")
	if answer.lower().startswith("y"):
		#Run lima
		print("***Run lima...***")
		time.sleep(2)
		ccs = args.ccs
		primers = args.primers
		
		lima = "lima --version"
		os.system(lima)
		
		lima_1 = "lima --same -m 75 -y 1 --dump-clips" + " " + ccs + " " + primers + " " + "output.same.bam"
		os.system(lima_1)
		lima_2 = "lima --diff -m 75 -y 1 --dump-clips" + " " + ccs + " " + primers + " " + "output.diff.bam"
		os.system(lima_2)
		
		#Run deconcat.py
		print("***Run deconcat.py...***")
		time.sleep(2)
		deconcat_1 = "python3 deconcat.py --method " + '"Jason-10X-3"' + " " + "output.same" + " " + "output.same.deconcat"
		os.system(deconcat_1)
		deconcat_2 = "python3 deconcat.py --method " + '"Jason-10X-3"' + " " + "output.diff" + " " + "output.diff.deconcat"
		os.system(deconcat_2)
		
		#Run extract_UMI_BC_post_deconcat.py
		print("***Run Run deconcat.py...***")
		time.sleep(2)
		extract_1 = "python3 extract_UMI_BC_post_deconcat.py -u 12 -b 16 --a3_clip_seq AGA --umi_type A3-10X output.same.deconcat.bam output.same.deconcat >output.same.deconcat.tagged.log"
		os.system(extract_1)
		extract_2 = "python extract_UMI_BC_post_deconcat.py -u 12 -b 16 --a3_clip_seq AGA --umi_type A3-10X output.diff.deconcat.bam output.diff.deconcat >output.diff.deconcat.tagged.log"
		os.system(extract_2)
		
		#Combine the same and diff outputs
		print("***Combine the same and diff outputs...***")
		time.sleep(2)
		combine = "python3 combine_same_diff.py -s output.same.deconcat.tagged -d output.diff.deconcat.tagged -o output.combined.flt"
		os.system(combine)
		
		subprocess.call("mv output.combined.flt flt.bam", shell = True)
		subprocess.call("ls flt.*", shell = True)
		time.sleep(2)
		
		#Run Refine
		print("***Run Refine...***")
		time.sleep(2)

		#Check the path of input file
		input_file2("flt.bam")

		#Current working directory
		current_path = os.path.abspath(os.getcwd())
		print(current_path)

		flt_bam = os.path.join(current_path, "flt.bam")
		primers = os.path.join(current_path, "primers.fasta")
		out = os.path.join(current_path, "fltnc.bam")
		refine = "isoseq3 refine" + " " + flt_bam + " " + primers + " " + out + " --require-polya"
		os.system(refine)

		subprocess.call("ls fltnc.*", shell = True)
		time.sleep(2)

		#Cluster Reads by Unique Founder Molecules
		print("***Cluster Reads by Unique Founder Molecules...***")
		time.sleep(2)
		
		#Check the path of input file
		input_file2("fltnc.bam")

		#Current working directory
		current_path = os.path.abspath(os.getcwd())
		print(current_path)
		
		fltnc_bam = os.path.join(current_path, "fltnc.bam")
		out = os.path.join(current_path, "dedup.bam")
		mismatch = ""
		mismatch = input("Enter maximum number of mismatches between tags: ")
		shift = ""
		shift = input("Enter tags may be shifted by at maximum of N bases: ")

		cluster = "isoseq3 dedup" + " " + fltnc_bam + " " + out + " --max-tag-mismatches " + str(mismatch) + " --max-tag-shift " + str(shift)
		os.system(cluster)
		
		subprocess.call("ls dedup.*", shell = True)
		time.sleep(2)

		#Make a dedup CSV file
		print("***Make a dedup CSV file...***")
		time.sleep(2)

		#Check the path of input file
		input_file2("dedup.fasta")

		make_csv = "python3 /opt/cDNA_Cupcake/cDNA_Cupcake/singlecell/make_csv_for_dedup.py"
		os.system(make_csv)
		subprocess.call("ls dedup.*", shell = True)
		time.sleep(2)

		#Align to Genome
		print("***Align to Genome...***")
		time.sleep(2)

		#Check the path of input file
		input_file2("dedup.fasta")

		#Add reference file
		reference = ""
		reference = input("Enter the path of the 'genome reference (.fa)': ")
		input_file2(reference)

		#Current working directory
		current_path = os.path.abspath(os.getcwd())
		print(current_path)

		dedup_fasta = os.path.join(current_path, "dedup.fasta")
		out = os.path.join(current_path, "dedup.fasta.sam")
		out2 = os.path.join(current_path, "dedup.fasta.sam.log")

		align = "minimap2 -t 30 -ax splice -uf --secondary=no -C5" + " " + reference + " " + dedup_fasta + " " + ">" + " " + out + " " + "2>" + " " + out2
		os.system(align)

		subprocess.call("ls dedup.*", shell = True)
		time.sleep(2)

		#Sort clustered.fasta.sam
		print("***Sort clustered.fasta.sam...***")
		sort = "sort -k 3,3 -k 4,4n" + " " + out + ">" + " " + "dedup.fasta.sorted.sam"
		os.system(sort)
		subprocess.call("ls dedup.*", shell = True)
		time.sleep(2)

		#Collapse into Unique Transcripts
		print("***Collapse into Unique Transcripts...***")
		time.sleep(2)

		#Check the path of input file
		input_file2("dedup.fasta")

		#Current working directory
		current_path = os.path.abspath(os.getcwd())
		print(current_path)

		dedup_fasta = os.path.join(current_path, "dedup.fasta")
		sam = os.path.join(current_path, "dedup.fasta.sorted.sam")
		out = os.path.join(current_path, "dedup.5merge")

		collapse = "collapse_isoforms_by_sam.py --input" + " " + dedup_fasta + " -s " + sam + " -c 0.99 -i 0.95 --gen_mol_count -o " + out
		os.system(collapse)

		subprocess.call("ls dedup.*", shell = True)
		time.sleep(2)

		#Compare Against Annotation
		exec(open("sqanti3.py").read())

		#Process into CSV Report and UMI/BC Error Correction
		print("***Process into CSV Report and UMI/BC Error Correction...***")
		time.sleep(2)

		#Current working directory
		current_path = os.path.abspath(os.getcwd())
		print(current_path)

		#Check file locations
		input_file2("dedup.5merge.collapsed.group.txt")
		input_file2("dedup.info.csv")
		input_file2(current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_classification.filtered_lite_classification.txt")

		#gzip dedup.info.csv
		gzip = "gzip dedup.info.csv"
		os.system(gzip)
		subprocess.call("ls *.gz", shell = True)
		time.sleep(2)

		group = current_path + "/" + "dedup.5merge.collapsed.group.txt"
		info = current_path + "/" + "dedup.info.csv.gz"
		lite_classification = current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_classification.filtered_lite_classification.txt"
		out = os.path.join(current_path, "dedup.annotated.csv")

		#Run collate_FLNC_gene_info.py
		gene_info = "python3 /opt/cDNA_Cupcake/cDNA_Cupcake/singlecell/collate_FLNC_gene_info.py" + " " + group + " " + info + " " + lite_classification + " " + out
		os.system(gene_info)
		subprocess.call("ls dedup.*", shell = True)
		time.sleep(2)

		#Producing a Seurat-compatible input from SQANTI3 isoforms
		print("***Producing a Seurat-compatible input from SQANTI3 isoforms...***")
		time.sleep(2)

		#Current working directory
		current_path = os.path.abspath(os.getcwd())
		print(current_path)

		#Check file locations
		input_file2("dedup.annotated.csv")
		input_file2(current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_classification.filtered_lite.gtf")

		annot = current_path + "/" + "dedup.annotated.csv"
		lite_classification = current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_classification.filtered_lite.gtf"
		out = os.path.join(current_path, "output")

		#Run make_seurat_input.py
		seurat = "python3 /opt/cDNA_Cupcake/cDNA_Cupcake/singlecell/make_seurat_input.py -i" + " " + annot + " -a " + lite_classification + " -o " + " " + out
		os.system(seurat)

		print(""" ***"<path_to>/output/isoforms_seurat" can be loaded into Seurat.***""")
		print("""***pbmc.data = Read10X(data.dir="<path_to>/output/isoforms_seurat").***""")
		time.sleep(2)

		#Linking molecule info to cell type, running IsoPhase, tagging BAMs
		exec(open("IsoPhase.py").read())
				
	elif answer.lower().startswith("n"):
		print("***Run extract_UMI_BC_post_deconcat.py...***")
		time.sleep(2)
		break



#Confirm primer sequences
primers_file = args.primers_file
print(primers_file.read())



#Run lima
print("***Run lima...***")
ccs = args.ccs
primers = args.primers
out = args.output

print("***Primer removal and demultiplexing***")
time.sleep(2)

lima = "lima --version"
os.system(lima)

lima = "lima --dump-clips" + " " + ccs + " " + primers + " " + out
os.system(lima)

subprocess.call("ls fl.*", shell = True)
time.sleep(2)



#Detect UMIs and Cell Barcodes
print("***Detect UMIs and Cell Barcodes...***")
time.sleep(2)

isoseq3 = "isoseq3 --version"
os.system(isoseq3)

#Change the input file name
subprocess.call("mv fl.bam lima.bam", shell = True)

#Check the path of input file
input_file2("lima.bam")

#Current working directory
current_path = os.path.abspath(os.getcwd())
print(current_path)

lima_bam = os.path.join(current_path, "lima.bam")
out = os.path.join(current_path, "flt.bam")

umi_cellBarcode = "isoseq3 tag" + " " + lima_bam + " " + out + " --design T-10U-16B"
os.system(umi_cellBarcode)

subprocess.call("ls flt.*", shell = True)
time.sleep(2)



#Run Refine
print("***Run Refine...***")
time.sleep(2)

#Check the path of input file
input_file2("flt.bam")

#Current working directory
current_path = os.path.abspath(os.getcwd())
print(current_path)

flt_bam = os.path.join(current_path, "flt.bam")
primers = os.path.join(current_path, "primers.fasta")
out = os.path.join(current_path, "fltnc.bam")
refine = "isoseq3 refine" + " " + flt_bam + " " + primers + " " + out + " --require-polya"
os.system(refine)

subprocess.call("ls fltnc.*", shell = True)
time.sleep(2)



#Cluster Reads by Unique Founder Molecules
print("***Cluster Reads by Unique Founder Molecules...***")
time.sleep(2)
		
#Check the path of input file
input_file2("fltnc.bam")

#Current working directory
current_path = os.path.abspath(os.getcwd())
print(current_path)
		
fltnc_bam = os.path.join(current_path, "fltnc.bam")
out = os.path.join(current_path, "dedup.bam")
mismatch = ""
mismatch = input("Enter maximum number of mismatches between tags: ")
shift = ""
shift = input("Enter tags may be shifted by at maximum of N bases: ")

cluster = "isoseq3 dedup" + " " + fltnc_bam + " " + out + " --max-tag-mismatches " + str(mismatch) + " --max-tag-shift " + str(shift)
os.system(cluster)
		
subprocess.call("ls dedup.*", shell = True)
time.sleep(2)



#Make a dedup CSV file
print("***Make a dedup CSV file...***")
time.sleep(2)

#Check the path of input file
input_file2("dedup.fasta")

make_csv = "python3 /opt/cDNA_Cupcake/cDNA_Cupcake/singlecell/make_csv_for_dedup.py"
os.system(make_csv)
subprocess.call("ls dedup.*", shell = True)
time.sleep(2)



#Align to Genome
print("***Align to Genome...***")
time.sleep(2)

#Check the path of input file
input_file2("dedup.fasta")

#Add reference file
reference = ""
reference = input("Enter the path of the 'genome reference (.fa)': ")
input_file2(reference)

#Current working directory
current_path = os.path.abspath(os.getcwd())
print(current_path)

dedup_fasta = os.path.join(current_path, "dedup.fasta")
out = os.path.join(current_path, "dedup.fasta.sam")
out2 = os.path.join(current_path, "dedup.fasta.sam.log")

align = "minimap2 -t 30 -ax splice -uf --secondary=no -C5" + " " + reference + " " + dedup_fasta + " " + ">" + " " + out + " " + "2>" + " " + out2
os.system(align)

subprocess.call("ls dedup.*", shell = True)
time.sleep(2)



#Sort clustered.fasta.sam
print("***Sort clustered.fasta.sam...***")
sort = "sort -k 3,3 -k 4,4n" + " " + out + ">" + " " + "dedup.fasta.sorted.sam"
os.system(sort)
subprocess.call("ls dedup.*", shell = True)
time.sleep(2)



#Collapse into Unique Transcripts
print("***Collapse into Unique Transcripts...***")
time.sleep(2)

#Check the path of input file
input_file2("dedup.fasta")

#Current working directory
current_path = os.path.abspath(os.getcwd())
print(current_path)

dedup_fasta = os.path.join(current_path, "dedup.fasta")
sam = os.path.join(current_path, "dedup.fasta.sorted.sam")
out = os.path.join(current_path, "dedup.5merge")

collapse = "collapse_isoforms_by_sam.py --input" + " " + dedup_fasta + " -s " + sam + " -c 0.99 -i 0.95 --gen_mol_count -o " + out
os.system(collapse)

subprocess.call("ls dedup.*", shell = True)
time.sleep(2)



#Compare Against Annotation
print("***Compare Against Annotation using SQANTI3...***")
time.sleep(2)

exec(open("sqanti3.py").read())



#Process into CSV Report and UMI/BC Error Correction
print("***Process into CSV Report and UMI/BC Error Correction...***")
time.sleep(2)

#Current working directory
current_path = os.path.abspath(os.getcwd())
print(current_path)

#Check file locations
input_file2("dedup.5merge.collapsed.group.txt")
input_file2("dedup.info.csv")
input_file2(current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_classification.filtered_lite_classification.txt")

#gzip dedup.info.csv
gzip = "gzip dedup.info.csv"
os.system(gzip)
subprocess.call("ls *.gz", shell = True)
time.sleep(2)

group = current_path + "/" + "dedup.5merge.collapsed.group.txt"
info = current_path + "/" + "dedup.info.csv.gz"
lite_classification = current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_classification.filtered_lite_classification.txt"
out = os.path.join(current_path, "dedup.annotated.csv")

#Run collate_FLNC_gene_info.py
gene_info = "python3 /opt/cDNA_Cupcake/cDNA_Cupcake/singlecell/collate_FLNC_gene_info.py" + " " + group + " " + info + " " + lite_classification + " " + out
os.system(gene_info)
subprocess.call("ls dedup.*", shell = True)
time.sleep(2)



#Producing a Seurat-compatible input from SQANTI3 isoforms
print("***Producing a Seurat-compatible input from SQANTI3 isoforms...***")
time.sleep(2)

#Current working directory
current_path = os.path.abspath(os.getcwd())
print(current_path)

#Check file locations
input_file2("dedup.annotated.csv")
input_file2(current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_classification.filtered_lite.gtf")

annot = current_path + "/" + "dedup.annotated.csv"
lite_classification = current_path + "/" + "sqanti3" + "/" + "dedup.5merge.collapsed_classification.filtered_lite.gtf"
out = os.path.join(current_path, "output")

#Run make_seurat_input.py
seurat = "python3 /opt/cDNA_Cupcake/cDNA_Cupcake/singlecell/make_seurat_input.py -i" + " " + annot + " -a " + lite_classification + " -o " + " " + out
os.system(seurat)

print(""" ***"<path_to>/output/isoforms_seurat" can be loaded into Seurat.***""")
print("""***pbmc.data = Read10X(data.dir="<path_to>/output/isoforms_seurat").***""")
time.sleep(2)



#Linking molecule info to cell type, running IsoPhase, tagging BAMs
exec(open("IsoPhase.py").read())