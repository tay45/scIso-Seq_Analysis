#!/usr/bin/python

import os
import subprocess
import sys
import time
import logging
import yaml

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Function to run subprocess commands with error checking
def run_command(command):
    try:
        logging.info(f"Running command: {command}")
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {e}")
        sys.exit(1)

# Load configuration from YAML file
def load_config(config_file):
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

# Main SQANTI3 pipeline
def run_sqanti3_pipeline(input_gtf, abundance, config, out_name="out_sqanti", output_dir="sqanti3"):
    logging.info("*** Running SQANTI3 Pipeline ***")
    time.sleep(2)

    if not os.path.exists(output_dir):
        logging.info(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir)

    # Build base command
    sqanti3_qc = f"sqanti3_qc.py {input_gtf} {config['ref_gtf']} {config['ref_genome']} -o {config['out_name']} -d {config['output_dir']} -fl {abundance} --cpus {config['cpus']} -n {config['chunks']} --genename --isoAnnotLite --gff3 {config['tappAS']} --report both"

    # Add optional flags only if provided in the config
    if config.get('cage_peak'):
        sqanti3_qc += f" --CAGE_peak {config['cage_peak']}"
    
    if config.get('polya_peak'):
        sqanti3_qc += f" --polyA_peak {config['polya_peak']}"
    
    if config.get('polya_motif'):
        sqanti3_qc += f" --polyA_motif_list {config['polya_motif']}"
    
    if config.get('coverage'):
        sqanti3_qc += f" -c {config['coverage']}"
    
    if config.get('short_reads'):
        sqanti3_qc += f" --short_reads {config['short_reads']}"

    run_command(sqanti3_qc)

    # Now run the SQANTI3 ML filter if specified
    if config.get('run_ml_filter', False):
        run_sqanti3_ml_filter(config)

    logging.info("*** Process Complete ***")

# Function to run SQANTI3 ML filter
def run_sqanti3_ml_filter(config):
    logging.info("*** Running SQANTI3 ML Filter ***")

    # Build the SQANTI3 ML filter command
    sqanti3_ml = f"sqanti3_filter.py ml {config['output_dir']}/{config['out_name']}_classification.txt --isoAnnotGFF3 {config['output_dir']}/{config['out_name']}.gff3 --gtf {config['output_dir']}/{config['out_name']}_corrected.gtf --faa {config['output_dir']}/{config['out_name']}_corrected.faa -o {config['out_name']}_ML -d {config['output_dir']}"

    # Add the optional --sam flag only if provided
    if config.get('sam'):
        sqanti3_ml += f" --sam {config['sam']}"

    # Add optional flags for the ML filter if provided in the config
    if config.get('percent_training'):
        sqanti3_ml += f" -t {config['percent_training']}"
    
    if config.get('tp'):
        sqanti3_ml += f" -p {config['tp']}"
    
    if config.get('tn'):
        sqanti3_ml += f" -n {config['tn']}"
    
    if config.get('threshold'):
        sqanti3_ml += f" -j {config['threshold']}"
    
    # Ensure that --force_fsm_in is properly passed with a boolean argument
    if config.get('force_fsm_in', False):
        sqanti3_ml += f" -f {config['force_fsm_in']}"

    if config.get('remove_columns'):
        sqanti3_ml += f" -r {config['remove_columns']}"

    if config.get('intrapriming'):
        sqanti3_ml += f" -i {config['intrapriming']}"

    if config.get('filter_mono_exonic', False):
        sqanti3_ml += " -e"

    # Add new options
    if config.get('max_class_size'):
        sqanti3_ml += f" -z {config['max_class_size']}"

    if config.get('intermediate_files', False):
        sqanti3_ml += " --intermediate_files"

    if config.get('version', False):
        sqanti3_ml += " -v"

    run_command(sqanti3_ml)
    logging.info("*** SQANTI3 ML Filter Completed ***")

# Function to run SQANTI3 rescue
def run_sqanti3_rescue(config):
    logging.info("*** Running SQANTI3 Rescue ***")
    
    # Build the SQANTI3 rescue command
    sqanti3_rescue = f"sqanti3_rescue.py ml {config['output_dir']}/{config['out_name']}_ML_MLresult_classification.txt --isoforms {config['output_dir']}/{config['out_name']}_corrected.fasta --gtf {config['output_dir']}/{config['out_name']}_ML.filtered.gtf -g {config['ref_gtf']} -f {config['ref_genome']} -k {config['ref_classif']} -o {config['out_name']}_rescue -d {config['output_dir']}"

    # Add optional flags for the rescue process
    if config.get('rescue_mono_exonic'):
        sqanti3_rescue += f" -e {config['rescue_mono_exonic']}"

    if config.get('mode'):
        sqanti3_rescue += f" --mode {config['mode']}"
    
    if config.get('randomforest'):
        sqanti3_rescue += f" -r {config['randomforest']}"
    
    if config.get('threshold'):
        sqanti3_rescue += f" -j {config['threshold']}"
    
    run_command(sqanti3_rescue)
    logging.info("*** SQANTI3 Rescue Completed ***")

# Main function to trigger SQANTI3
if __name__ == "__main__":
    # Load configuration file
    config = load_config("config.yaml")

    # Initial files
    input_gtf = config.get("input_gtf", "collapsed.gtf")
    abundance = config.get("abundance", "collapsed.abundance.tsv")

    # Check if there are chained jobs (before running SQANTI3-QC)
    chain_jobs_value = str(config.get('chain_jobs', 'no')).lower()
    
    if chain_jobs_value == 'yes':
        logging.info("*** Chaining jobs enabled. Running chain.py script ***")
        # Run the chain.py script with the same config file to process chaining
        chain_command = f"python3 chain.py config.yaml"
        run_command(chain_command)

        # Check if chaining outputs exist and update input files accordingly
        chained_gtf = "all_samples.chained.gff"
        chained_abundance = "all_samples.chained_count.txt"

        if os.path.exists(chained_gtf) and os.path.exists(chained_abundance):
            logging.info(f"Chaining completed, using chained files: {chained_gtf} and {chained_abundance}")
            input_gtf = chained_gtf
            abundance = chained_abundance
        else:
            logging.error("Chaining output files not found! Ensure chaining ran successfully.")
            sys.exit(1)
    else:
        logging.info("*** Skipping chaining as per user's choice ***")

    # Run the SQANTI3 pipeline
    run_sqanti3_pipeline(input_gtf, abundance, config)

    # Run SQANTI3 Rescue if specified
    if config.get('run_rescue', False):
        run_sqanti3_rescue(config)

    logging.info("*** All processes completed ***")
