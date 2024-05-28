#SNP (single nucleotide polymorphism)
#INDEL (insertion/deletion)
#MIXED (combination of SNPs and indels at a single position)
#MNP (multi-nucleotide polymorphism, e.g. a dinucleotide substitution)

import os
import subprocess
import pandas as pd
import sys

def generate_bcftools_stats_for_folder(vcf_folder, output_folder):
    """Generate statistics from VCF files in a folder using bcftools stats."""
    os.makedirs(output_folder, exist_ok=True)
    for vcf_file in os.listdir(vcf_folder):
        if vcf_file.endswith(".vcf"):
            sample_name = vcf_file.split("-")[0]
            vcf_file_path = os.path.join(vcf_folder, vcf_file)
            output_file = os.path.join(output_folder, f"{sample_name}_bcftools_stats.txt")
            generate_bcftools_stats(vcf_file_path, output_file)

def generate_bcftools_stats(vcf_file_path, output_file):
    """Generate statistics from a VCF file using bcftools stats."""
    command = f"bcftools stats {vcf_file_path}"
    try:
        result = subprocess.run(command, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        with open(output_file, "w") as f:
            f.write(result.stdout)
    except Exception as e:
        print(f"An error occurred: {str(e)}")

def parse_stats_files(output_folder):
    """Parse statistics files and create a summary TSV file."""
    headers = ["number of records", "number of SNPs", "number of MNPs", "number of indels", "number of others"]
    summary_data = {header: [] for header in headers}
    filenames = []

    for stat_file in os.listdir(output_folder):
        if stat_file.endswith("_bcftools_stats.txt"):
            sample_name = stat_file.split("_")[0]
            filenames.append(sample_name)
            file_path = os.path.join(output_folder, stat_file)
            data = parse_specific_stats(file_path)
            for header in headers:
                summary_data[header].append(data.get(header, '0'))  # Default to '0' if not found

    # Create DataFrame
    df = pd.DataFrame(summary_data, index=filenames)
    # Transpose DataFrame to meet the expected format
    df = df.T

    # Rename row names
    new_row_names = ["number of variants", "number of SNPs", "number of MNPs", "number of indels", "number of complex"]
    df.index = new_row_names

    # Write to TSV
    df.to_csv(os.path.join(output_folder, "summary_stats.tsv"), sep='\t')

def parse_specific_stats(file_path):
    metrics = {
        "number of records": "SN\t0\tnumber of records:",
        "number of SNPs": "SN\t0\tnumber of SNPs:",
        "number of MNPs": "SN\t0\tnumber of MNPs:",
        "number of indels": "SN\t0\tnumber of indels:",
        "number of others": "SN\t0\tnumber of others:"
    }
    results = {}
    with open(file_path, "r") as file:
        for line in file:
            for key, value in metrics.items():
                if value in line:
                    results[key] = line.split()[-1]
    return results



vcf_folder = sys.argv[1]
output_folder = sys.argv[2]
generate_bcftools_stats_for_folder(vcf_folder, output_folder)
parse_stats_files(output_folder)
