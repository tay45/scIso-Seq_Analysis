scIso-Seq Analysis Pipeline

This single cell (sc) Iso-Seq pipeline is based on the PacBio IsoSeq3 v3.4.0 (https://github.com/PacificBiosciences/IsoSeq), SQANTI3 v4.2 (https://github.com/ConesaLab/SQANTI3) and the tools of the cDNA_Cupcake v28.0.0 (https://github.com/Magdoll/cDNA_Cupcake). It can help the people, who don't have the prior experience of the PacBio scIso-Seq data analysis, to easily proceed that using a single command line, and the interactive input function in the pipeline can also guide them to choose the right next step of the analysis.

It is designed for executing in the Apollo server (apollo-acc.coh.org) of the City of Hope (Duarte, CA) using the interactive job mode (http://apollo.coh.org/user-guide/interactivejobs/). Before executing the analysis, the below modules must be loaded in advance to the server.

    module load smrtlink/11.0
    module load cDNA_Cupcake
    module load SQANTI3
    module load gffread/0.12.7

To execute the pipeline, you need to download the python scripts of the pipeline (git clone https://github.com/tay45/scIso-Seq_Analysis.git), first. And, copy the scripts to the folder to run the analysis (../running_folder/scIsoSeq3.py, sqanti3.py and IsoPhase.py). The folder should contain raw read (ccs.bam), index (ccs.bam.pbi) and primer (primers.fasta) files. Finally, download the IsoAnnot script (https://isoannot.tappas.org/isoannot-lite/), IsoAnnotLite_v2.7.0_SQ3.py, and copy it to the same folder.

And, execute the below command line in the shell.

    python3 scIsoSeq3.py -c ccs.bam -p primers.fasta -o fl.bam -pf primers.fasta

The input must be HiFi reads (ccs.bam; QV>20). If the raw data is the continuous long reads (CLRs) in movieX.subreads.bam, please produce the ccs.bam via the below command line to generate the correct input.

    ccs movieX.subreads.bam movieX.ccs.bam --min-rq 0.9

During the running of the pipeline, it asks your intention of the analysis (e.g., whether or not you want to provide short-read fastq file for SQANTI annotation) to decide the direction the next step or ask you to enter a certain file (e.g., a genome reference) to proceed the next analysis.

The 'supporting_files' folder (/net/isi-dcnl/ifs/user_data/Seq/PacBio/thkang/Pipelines/IsoSeq3_Sqanti3/supporting_files) contains the files requiring for the SQANTI3 execution. If you need a more information regarding the files, please refer to the instruction of the SQANTI3 (https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-Quality-Control).

    Genome reference: GRCh38.p12.genome.fa; hg38.fa
    Annotation: gencode.v30.chr_patch_hapl_scaff.annotation.gtf; gencode.v30.annotation.gtf
    CAGE_Peak: hg38.cage_peak_phase1and2combined_coord.bed
    polyA_Peak: atlas.clusters.2.0.GRCh38.96.bed
    polyA_Motif: mouse_and_human.polyA_motif.txt
    tappAS-annotation: Homo_sapiens_GRCh38_Ensembl_86.gff3; Homo_sapiens_GRCh38_RefSeq_78.gff3
    Intropolis Junction BED: intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified

Regarding the outputs;

    sqanti3 folder: containing the output of the SQANTI3 analysis using the transcrpts from a single sample.
    colored.CPM.bed12: making colored BED format shaded by isoform abundance.
    isoannot.gff3: An input for the tappAS analysis.
    dedup.annotated.csv: generated a collated CSV file that links each mapped FLNC read to its classified genes and transcripts based on SQANTI3's output. 
    output folder: generated a Seurat-compatible input (pbmc.data = Read10X(data.dir="<path_to>/output/isoforms_seurat").
    dedup.annotated.celltype.csv: a file added the celltype column to the annotated CSV file.
    by_loci folder: containing information about the gene loci with sufficient coverage for phasing (by_loci/PB.X_sizeN).
    by_loci/PB.X_sizeN/phased.nopartial.cleaned.vcf: phasing using reads that cover every SNP base.
    by_loci/PB.X_sizeN/phased.partial.cleaned.vcf: phasing using reads that cover at least one SNP base.
    by_loci/PB.X_sizeN/ccs.hg38.sorted.tagged.bam: the aligned read BAM file with phasing information in the RG tag which can be visualized by IGV.
    summarized.isophase_results.txt: a summary file which shows the 'num_snp', the 'num_hap_nopartial' and all 'num_hap_withpartial'.
    IsoSeq_IsoPhase.vcf: a summarized SNP VCF file collecting all the VCFs from each separate directory.

For executing this pipeline in other platforms, please satisfy the prerequisite through installing below tools.

    IsoSeq3 (https://github.com/PacificBiosciences/IsoSeq)
    SQANTI3 and cDNA_Cupcake (https://github.com/ConesaLab/SQANTI3/wiki/SQANTI3-dependencies-and-installation)
    gffread (https://anaconda.org/bioconda/gffread)

