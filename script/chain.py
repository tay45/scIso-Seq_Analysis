#!/usr/bin/python

import os
import sys
import yaml
import glob

# Function to check for files by pattern
def find_file_by_pattern(directory, pattern):
    files = glob.glob(os.path.join(directory, pattern))
    if files:
        return os.path.basename(files[0])  # Return only the file name
    else:
        print(f"File matching pattern '{pattern}' not found in {directory}. Please check the file location.")
        sys.exit(1)

# Load the configuration from the config.yaml file
def load_config(config_file):
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

# Generate the sample.config file based on config.yaml
def create_sample_config(samples):
    sample_config_path = "sample.config"
    
    with open(sample_config_path, 'w') as config:
        for sample in samples:
            # Find the required files by pattern in each sample directory
            collapsed_group = find_file_by_pattern(sample['path'], '*.group.txt')
            collapsed_gff = find_file_by_pattern(sample['path'], '*.gff')
            collapsed_abundance = find_file_by_pattern(sample['path'], '*.abundance.txt')
            
            # Write the SAMPLE lines with the sample name and path
            config.write(f"SAMPLE={sample['name']};{sample['path']}\n")
        
        # Write the filenames (without paths) to the config file
        config.write(f"GROUP_FILENAME={collapsed_group}\n")
        config.write(f"GFF_FILENAME={collapsed_gff}\n")
        config.write(f"COUNT_FILENAME={collapsed_abundance}\n")
        
        # Optionally include FASTQ_FILENAME if a fastq file exists in the directory
        fastq_file = glob.glob(os.path.join(sample['path'], '*.fastq'))
        if fastq_file:
            config.write(f"FASTQ_FILENAME={os.path.basename(fastq_file[0])}\n")

    print(f"Configuration file created: {sample_config_path}")
    return sample_config_path

# Main function to generate the config file and run chain_samples.py
def main():
    config = load_config("config.yaml")
    samples = config.get("samples", [])

    if not samples:
        print("No samples found in configuration. Please specify samples under 'samples' in config.yaml.")
        sys.exit(1)

    # Create the sample config file
    sample_config = create_sample_config(samples)
    
    # Run the chain_samples.py script with the generated sample.config
    chain_command = f"chain_samples.py {sample_config} count_fl --dun-merge-5-shorter"
    print(f"Executing command: {chain_command}")
    os.system(chain_command)

if __name__ == "__main__":
    main()
