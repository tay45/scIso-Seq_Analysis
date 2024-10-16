import argparse
import subprocess
import sys
import logging
import os
import resource

# Set up logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(message)s',
                    handlers=[
                        logging.FileHandler("isoseq_pipeline.log"),
                        logging.StreamHandler(sys.stdout)
                    ])

# Run a command with logging for stdout and stderr
def run_command(command, description):
    logging.info(f"*** {description} ***")
    logging.info(f"Running command: {command}")
    
    try:
        # Use universal_newlines instead of text for Python < 3.7 compatibility
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

        # Stream output from the command in real-time
        for stdout_line in iter(process.stdout.readline, ''):
            logging.info(stdout_line.strip())
        for stderr_line in iter(process.stderr.readline, ''):
            logging.error(stderr_line.strip())

        process.stdout.close()
        process.stderr.close()
        returncode = process.wait()

        if returncode != 0:
            logging.error(f"Command failed with exit code {returncode}")
            raise subprocess.CalledProcessError(returncode, command)
        logging.info(f"Finished {description} successfully")

    except subprocess.CalledProcessError as e:
        logging.error(f"Command {command} failed: {e}")
        sys.exit(1)

# Set explicit resource limits
def set_resource_limits():
    # Setting memory limit to avoid excessive memory usage crashes (customize as per your environment)
    try:
        soft, hard = resource.getrlimit(resource.RLIMIT_AS)
        resource.setrlimit(resource.RLIMIT_AS, (soft, hard))
        logging.info(f"Memory resource limit set to: {soft}")
    except ValueError as e:
        logging.error(f"Error setting resource limits: {e}")

# Check if file exists
def check_file_exists(filepath):
    if not os.path.isfile(filepath):
        logging.error(f"File not found: {filepath}")
        raise FileNotFoundError(f"File not found: {filepath}")

# Define the function for skera split
def skera_split(consensusreadset, fragment_adapters, num_threads):
    logging.info("Starting skera split")
    
    output_bam = "segmented.consensusreadset.bam"  # Fixed output file
    logging.info(f"Output BAM file will be: {output_bam}")
    
    skera_command = f"skera split --log-level INFO --log-file skera.log --alarms alarms.json -j {num_threads} \"{consensusreadset}\" \"{fragment_adapters}\" \"{output_bam}\""
    
    run_command(skera_command, f"Running skera split with output to {output_bam}")
    
    return output_bam

# Define the function for lima
def lima(segmented_bam, cdna_barcodes, num_threads):
    logging.info("Starting lima")
    
    scisoseq_bam = "scisoseq.bam"  # Fixed output file name
    logging.info(f"Output BAM file will be: {scisoseq_bam}")
    
    lima_command = f"lima --log-level DEBUG --log-file lima_isoseq.log -j {num_threads} --alarms alarms.json --isoseq --per-read \"{segmented_bam}\" \"{cdna_barcodes}\" \"{scisoseq_bam}\""
    
    run_command(lima_command, f"Running lima with output to {scisoseq_bam}")
    
    return scisoseq_bam

# Define the function for isoseq tag
def isoseq_tag(scisoseq_bam, design, num_threads):
    logging.info("Starting isoseq tag")
    
    tagged_bam = "scisoseq.5p--3p.tagged.bam"
    tag_command = f"isoseq tag --log-level INFO --log-file isoseq_tag.log -j {num_threads} --alarms alarms.json --design {design} \"scisoseq.5p--3p.bam\" \"{tagged_bam}\""
    
    run_command(tag_command, f"Running isoseq tag with output to {tagged_bam}")
    
    return tagged_bam

# Define the function for isoseq refine
def isoseq_refine(tagged_bam, cdna_barcodes, num_threads):
    logging.info("Starting isoseq refine")
    
    refined_bam = "scisoseq.5p--3p.tagged.refined.bam"
    refine_command = f"isoseq refine --log-level INFO --log-file isoseq_refine.log -j {num_threads} --alarms alarms.json --require-polya \"{tagged_bam}\" \"{cdna_barcodes}\" \"{refined_bam}\""
    
    run_command(refine_command, f"Running isoseq refine with output to {refined_bam}")
    
    return refined_bam

# Define the function for isoseq correct
def isoseq_correct(refined_bam, barcodes_gz, num_threads):
    logging.info("Starting isoseq correct")
    
    corrected_bam = "scisoseq.5p--3p.tagged.refined.corrected.bam"
    correct_command = f"isoseq correct --log-level DEBUG --verbose --log-file isoseq_correct.log -j {num_threads} --alarms alarms.json --method knee --percentile 99 --barcodes \"{barcodes_gz}\" \"{refined_bam}\" \"{corrected_bam}\""
    
    run_command(correct_command, f"Running isoseq correct with output to {corrected_bam}")
    
    return corrected_bam

# Define the function for samtools sort
def samtools_sort(corrected_bam, num_threads):
    logging.info("Starting samtools sort")
    
    sorted_bam = "scisoseq.5p--3p.tagged.refined.corrected.sorted.bam"
    sort_command = f"samtools sort -t CB -@ {num_threads} -m 768M \"{corrected_bam}\" -o \"{sorted_bam}\""
    
    run_command(sort_command, f"Running samtools sort with output to {sorted_bam}\"")
    
    return sorted_bam

# Define the function for isoseq groupdedup
def isoseq_groupdedup(sorted_bam, num_threads):
    logging.info("Starting isoseq groupdedup")
    
    dedup_bam = "scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.bam"
    groupdedup_command = f"isoseq groupdedup --log-level INFO --log-file isoseq_dedup.log --alarms alarms.json -j {num_threads} \"{sorted_bam}\" \"{dedup_bam}\""
    
    run_command(groupdedup_command, f"Running isoseq groupdedup with output to {dedup_bam}\"")
    
    return dedup_bam

# Define the function for isoseq bcstats
def isoseq_bcstats(sorted_bam, num_threads):
    logging.info("Starting isoseq bcstats")
    
    bcstats_command = f"isoseq bcstats --log-level INFO --log-file isoseq_bcstats.log -j {num_threads} --alarms alarms.json --method knee --percentile 99 -o bcstats_report.tsv --json bcstats.report.json \"{sorted_bam}\""
    
    run_command(bcstats_command, f"Running isoseq bcstats on {sorted_bam}")
    
    return "bcstats_report.tsv"

# Define the function for pbmm2 align
def pbmm2_align(reference_fasta, dedup_bam, num_threads):
    logging.info("Starting pbmm2 align")
    
    mapped_bam = "scisoseq.mapped.bam"
    pbmm2_command = f"pbmm2 align --log-level INFO --log-file pbmm2.log -j {num_threads} --alarms alarms.json --preset ISOSEQ --sort --report-json mapping_stats.report.json \"{reference_fasta}\" \"{dedup_bam}\" \"{mapped_bam}\""
    
    run_command(pbmm2_command, f"Running pbmm2 align with output to {mapped_bam}")
    
    return mapped_bam

# Define the function for isoseq collapse
def isoseq_collapse(mapped_bam, num_threads):
    logging.info("Starting isoseq collapse")
    
    collapsed_gff = "scisoseq.mapped_transcripts.collapse.gff"
    collapse_command = f"isoseq collapse --log-level INFO --log-file isoseq_collapse.log -j {num_threads} --alarms alarms.json \"{mapped_bam}\" \"{collapsed_gff}\""
    
    run_command(collapse_command, f"Running isoseq collapse with output to {collapsed_gff}")
    
    return collapsed_gff

# Define the function for pigeon_prepare
def pigeon_prepare(collapsed_gff, annotation_gtf):
    logging.info("Starting pigeon_prepare")
    
    pigeon_command = f"pigeon prepare --log-level INFO --log-file pigeon_prepare.log \"{collapsed_gff}\" \"{annotation_gtf}\""
    
    run_command(pigeon_command, f"Running pigeon_prepare with inputs {collapsed_gff}, {annotation_gtf}")
    
    return "pigeon_prepare completed"

# Define the function for pigeon classify
def pigeon_classify(flnc_abundance, sorted_gff, annotation_gtf, reference_fasta):
    logging.info("Starting pigeon classify")
    
    # Modify the annotation GTF filename to add '.sorted'
    sorted_gtf = annotation_gtf.replace(".gtf", ".sorted.gtf")
    
    # Check file existence
    check_file_exists(sorted_gff)
    check_file_exists(flnc_abundance)
    check_file_exists(reference_fasta)
    check_file_exists(sorted_gtf)  # Ensure the sorted GTF file exists

    # Run classify command
    classify_command = f"pigeon classify --log-level INFO --log-file pigeon-classify.log --out-dir ./ --out-prefix scisoseq_pigeon_classify --flnc {flnc_abundance} {sorted_gff} {sorted_gtf} {reference_fasta}"
    
    run_command(classify_command, f"Running pigeon classify with inputs {flnc_abundance}, {sorted_gff}, {sorted_gtf}, {reference_fasta}")

# Define the function for pigeon filter
def pigeon_filter(classification_txt, isoforms_gff):
    logging.info("Starting pigeon filter")

    # Correct order of arguments in the filter command
    filter_command = f"pigeon filter {classification_txt} --isoforms {isoforms_gff} --log-level INFO --log-file pigeon-filter.log"
    
    run_command(filter_command, f"Running pigeon filter with inputs {classification_txt}, {isoforms_gff}")

    return "pigeon filter completed"

# Define the function for pigeon make-seurat
def pigeon_make_seurat(annotation_gtf, dedup_fasta, group_txt, classification_txt, num_threads):
    logging.info("Starting pigeon make-seurat")
    
    # Modify the annotation GTF filename to add '.sorted'
    sorted_gtf = annotation_gtf.replace(".gtf", ".sorted.gtf")
    
    # Check file existence
    check_file_exists(sorted_gtf)
    check_file_exists(dedup_fasta)
    check_file_exists(group_txt)
    check_file_exists(classification_txt)

    # Run make-seurat command
    make_seurat_command = f"pigeon make-seurat --log-level INFO --log-file pigeon-make-seurat.log --num-threads {num_threads} --annotations {sorted_gtf} --dedup {dedup_fasta} --group {group_txt} --out-dir . --out-prefix scisoseq_pigeon_seurat {classification_txt}"
    
    run_command(make_seurat_command, f"Running pigeon make-seurat with inputs {sorted_gtf}, {dedup_fasta}, {group_txt}, {classification_txt}")

# Main function for argument parsing and execution
def main():
    parser = argparse.ArgumentParser(description="Run the Iso-Seq pipeline")

    # Define the command-line arguments
    parser.add_argument("-c", "--consensusreadset", required=True, help="Path to the consensus read set XML file")
    parser.add_argument("-fa", "--fragment_adapters", required=False, help="Path to the fragment adapters FASTA file")  # Now optional
    parser.add_argument("-bc", "--cdna_barcodes", required=True, help="Path to the cDNA barcodes FASTA file for refine and correct")
    parser.add_argument("-bcgz", "--barcodes_gz", required=True, help="Path to the gzipped barcodes text file")
    parser.add_argument("-d", "--design", required=True, help="Design for isoseq tag (e.g., 16B-10U-10X-T)")
    parser.add_argument("-ref", "--reference_fasta", required=True, help="Path to the reference FASTA file")
    parser.add_argument("-gtf", "--annotation_gtf", required=True, help="Path to the annotation GTF file")
    parser.add_argument("-j", "--num_threads", required=False, default=7, type=int, help="Number of threads to use (default is 7)")

    # Parse arguments
    args = parser.parse_args()

    # Set resource limits
    set_resource_limits()

    # Only run skera split if the -fa (fragment adapters) flag is provided
    if args.fragment_adapters:
        segmented_bam = skera_split(args.consensusreadset, args.fragment_adapters, args.num_threads)
    else:
        # If -fa is not provided, assume the input BAM is already segmented
        segmented_bam = args.consensusreadset

    # Run lima
    scisoseq_bam = lima(segmented_bam, args.cdna_barcodes, args.num_threads)

    # Run isoseq tag
    tagged_bam = isoseq_tag(scisoseq_bam, args.design, args.num_threads)

    # Run isoseq refine
    refined_bam = isoseq_refine(tagged_bam, args.cdna_barcodes, args.num_threads)

    # Run isoseq correct
    corrected_bam = isoseq_correct(refined_bam, args.barcodes_gz, args.num_threads)

    # Run samtools sort
    sorted_bam = samtools_sort(corrected_bam, args.num_threads)

    # Run isoseq groupdedup
    dedup_bam = isoseq_groupdedup(sorted_bam, args.num_threads)

    # Run isoseq bcstats
    bcstats_report = isoseq_bcstats(sorted_bam, args.num_threads)

    # Run pbmm2 align
    mapped_bam = pbmm2_align(args.reference_fasta, dedup_bam, args.num_threads)

    # Run isoseq collapse
    collapsed_gff = isoseq_collapse(mapped_bam, args.num_threads)

    # Run pigeon_prepare
    pigeon_prepare_result = pigeon_prepare(collapsed_gff, args.annotation_gtf)

    # Run pigeon classify
    pigeon_classify_result = pigeon_classify("scisoseq.mapped_transcripts.collapse.abundance.txt", "scisoseq.mapped_transcripts.collapse.sorted.gff", args.annotation_gtf, args.reference_fasta)

    # Run pigeon filter
    pigeon_filter_result = pigeon_filter("scisoseq_pigeon_classify_classification.txt", "scisoseq.mapped_transcripts.collapse.sorted.gff")

    # Run pigeon make-seurat
    pigeon_make_seurat(args.annotation_gtf, "scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta", "scisoseq.mapped_transcripts.collapse.group.txt", "scisoseq_pigeon_classify_classification.filtered_lite_classification.txt", args.num_threads)

    logging.info(f"Pipeline completed successfully, collapsed GFF: {collapsed_gff}, BCStats report: {bcstats_report}, Pigeon Prepare: {pigeon_prepare_result}, Pigeon Classify: {pigeon_classify_result}, Pigeon Filter: {pigeon_filter_result}, Pigeon Make-Seurat completed")

if __name__ == "__main__":
    main()