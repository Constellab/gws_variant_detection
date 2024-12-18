import os
import pandas as pd
import subprocess
import numpy as np
import re
import sys
import pysam
import shutil
from multiprocessing import Pool

class VCFWorkflow:
    def __init__(self, bam_folder, ref_genome, gene_name, gene_positions, minipileup_path, output_folder,
                 counts_matrices_folder, frequencies_matrices_folder, known_variants_file, threads):
        self.bam_folder = bam_folder
        self.ref_genome = ref_genome
        self.gene_name = gene_name
        self.gene_positions = gene_positions  # List of positions of interest
        self.minipileup_path = minipileup_path
        self.output_folder = output_folder
        self.counts_matrices_folder = counts_matrices_folder
        self.frequencies_matrices_folder = frequencies_matrices_folder
        self.known_variants_file = known_variants_file
        self.threads = threads  # Number of threads for parallel processing
        self.create_folder(self.output_folder)
        self.create_folder(self.counts_matrices_folder)
        self.create_folder(self.frequencies_matrices_folder)
        self.per_read_bam_folder = os.path.join(self.output_folder, "per_read_bams")
        self.per_read_vcf_folder = os.path.join(self.output_folder, "per_read_vcfs")
        self.normalized_vcf_folder = os.path.join(self.output_folder, "normalized_vcfs")
        self.variant_names = self.get_variant_names()
        self.barcode_list = self.get_barcode_list()
        self.total_reads_per_barcode = {}  # Store total reads per barcode

    def create_folder(self, folder_path):
        """Ensure a folder exists; if not, create it."""
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            print(f"Created folder: {folder_path}")

    def get_variant_names(self):
        """Get mapping from (CHROM, POS, REF, ALT) to variant names for the selected gene."""
        # Only include variants for the selected gene
        all_variant_names = {
            ('CYP51', 149, 'T', 'C'): 'L50S',
            ('CYP51', 320, 'A', 'T'): 'D107V',
            ('CYP51', 401, 'A', 'G'): 'D134G',
            ('CYP51', 407, 'T', 'C'): 'V136A',
            ('CYP51', 410, 'A', 'T'): 'Y137F',
            ('CYP51', 433, 'A', 'C'): 'M145L_CTG',
            ('CYP51', 433, 'A', 'T'): 'M145L_TTG',
            ('CYP51', 590, 'A', 'G'): 'N178S',
            ('CYP51', 620, 'G', 'A'): 'S188N',
            ('CYP51', 679, 'T', 'A'): 'S208T',
            ('CYP51', 841, 'C', 'T'): 'H262Y',
            ('CYP51', 1398, 'T', 'C'): 'C372R',
            ('CYP51', 1408, 'C', 'G'): 'A379G',
            ('CYP51', 1413, 'A', 'G'): 'I381V',
            ('CYP51', 1646, 'CTAT', 'C'): 'Y459',
            ('CYP51', 1648, 'A', 'C'): 'Y459S',
            ('CYP51', 1649, 'TGGC', 'T'): 'G460',
            ('CYP51', 1653, 'T', 'C'): 'Y461H',
            ('CYP51', 1654, 'A', 'C'): 'Y461S',
            ('CYP51', 1811, 'C', 'A'): 'N513K_AAA',
            ('CYP51', 1811, 'C', 'G'): 'N513K_AAG',
            ('CYP51', 1854, 'T', 'A'): 'S524T',
            ('CytB', 110, 'G', 'T'): 'G37V',
            ('CytB', 428, 'G', 'C'): 'G143A',
            ('SdhB', 788, 'A', 'C'): 'N225T',
            ('SdhB', 788, 'A', 'T'): 'N225I',
            ('SdhB', 908, 'G', 'C'): 'R265P',
            ('SdhB', 914, 'A', 'T'): 'H267L',
            ('SdhB', 917, 'C', 'T'): 'T268I',
            ('SdhB', 916, 'A', 'G'): 'T268A',
            ('SdhC', 365, 'C', 'A'): 'T79N',
            ('SdhC', 368, 'G', 'C'): 'W80S',
            ('SdhC', 379, 'G', 'T'): 'A84F',
            ('SdhC', 386, 'A', 'G'): 'N86S',
            ('SdhC', 387, 'C', 'A'): 'N86K_AAA',
            ('SdhC', 387, 'C', 'G'): 'N86K_AAG',
            ('SdhC', 508, 'C', 'G'): 'P127A',
            ('SdhC', 582, 'G', 'C'): 'R151S_AGC',
            ('SdhC', 582, 'G', 'T'): 'R151S_AGT',
            ('SdhC', 581, 'G', 'C'): 'R151T',
            ('SdhC', 581, 'G', 'T'): 'R151M',
            ('SdhC', 580, 'A', 'G'): 'R151G',
            ('SdhC', 584, 'A', 'G'): 'H152R',
            ('SdhC', 611, 'T', 'G'): 'I161S',
            ('SdhC', 625, 'G', 'A'): 'V166M',
            ('SdhC', 632, 'C', 'G'): 'T168R',
            ('SdhD', 139, 'C', 'T'): 'R47W',
            ('SdhD', 148, 'A', 'T'): 'I50F',
            ('SdhD', 148, 'A', 'C'): 'I50L',
            ('SdhD', 399, 'A', 'G'): 'M114V',
            ('SdhD', 497, 'A', 'G'): 'D129G',
            ('SdhD', 498, 'C', 'A'): 'D129E_GAA',
            ('SdhD', 498, 'C', 'G'): 'D129E_GAG',
            # Add more mappings for other genes if needed
        }
        variant_names = {k: v for k, v in all_variant_names.items() if k[0] == self.gene_name}
        return variant_names

    def get_barcode_list(self):
        """Get list of barcodes from BAM filenames."""
        barcodes = set()
        for bam_file in os.listdir(self.bam_folder):
            if bam_file.endswith('.bam'):
                barcode = bam_file.split('-')[0]  # Assuming barcode is before '-'
                barcodes.add(barcode)
        return list(barcodes)

    def parse_known_variants_file(self):
        """Parse the known_variants_file to get variant counts per barcode."""
        known_variants_df = pd.read_csv(self.known_variants_file, sep='\t')
        # Expected columns: 'variant_name', 'barcodeXX_VARIANT_COUNT', 'barcodeXX_VARIANT_FREQ', etc.
        variant_counts_per_barcode = {}
        # Extract all barcode columns ending with '_VARIANT_COUNT'
        barcode_count_columns = [col for col in known_variants_df.columns if col.endswith('_VARIANT_COUNT')]
        for barcode_col in barcode_count_columns:
            barcode = barcode_col.replace('_VARIANT_COUNT', '')
            if barcode in self.barcode_list:
                counts = known_variants_df[['variant_name', barcode_col]].set_index('variant_name')[barcode_col].to_dict()
                variant_counts_per_barcode[barcode] = counts
            else:
                variant_counts_per_barcode[barcode] = {}
        return variant_counts_per_barcode

    def parse_known_variants_file_frequencies(self):
        """Parse the known_variants_file to get variant frequencies per barcode."""
        known_variants_df = pd.read_csv(self.known_variants_file, sep='\t')
        variant_freqs_per_barcode = {}
        barcode_freq_columns = [col for col in known_variants_df.columns if col.endswith('_VARIANT_FREQ')]
        for barcode_col in barcode_freq_columns:
            barcode = barcode_col.replace('_VARIANT_FREQ', '')
            if barcode in self.barcode_list:
                freqs = known_variants_df[['variant_name', barcode_col]].set_index('variant_name')[barcode_col].to_dict()
                # Remove '%' if present and convert to float
                freqs = {k: float(str(v).replace('%', '').replace(',', '.')) for k, v in freqs.items()}
                variant_freqs_per_barcode[barcode] = freqs
            else:
                variant_freqs_per_barcode[barcode] = {}
        return variant_freqs_per_barcode

    def split_bam_by_read(self):
        """Split each BAM file into multiple BAM files, each containing a single read."""
        self.create_folder(self.per_read_bam_folder)
        bam_files = [f for f in os.listdir(self.bam_folder) if f.endswith('.bam')]
        bam_paths = [os.path.join(self.bam_folder, bam_file) for bam_file in bam_files]
        with Pool(processes=self.threads) as pool:
            pool.map(self.process_bam_file, bam_paths)

    def process_bam_file(self, bam_path):
        """Process a single BAM file to split it into per-read BAM files."""
        bam_file = os.path.basename(bam_path)
        bam_basename = os.path.splitext(bam_file)[0]
        barcode = bam_basename.split('-')[0]
        barcode_folder = os.path.join(self.per_read_bam_folder, barcode)
        self.create_folder(barcode_folder)
        print(f"Processing BAM file: {bam_file}")

        bamfile = pysam.AlignmentFile(bam_path, "rb")
        read_counter = 0

        for read in bamfile.fetch(until_eof=True):
            read_counter += 1
            per_read_bam_path = os.path.join(barcode_folder, f"{bam_basename}_read{read_counter}.bam")
            with pysam.AlignmentFile(per_read_bam_path, "wb", template=bamfile) as out_bam:
                out_bam.write(read)
            # Index the per-read BAM file
            pysam.index(per_read_bam_path)

        bamfile.close()
        print(f"Split {bam_file} into {read_counter} per-read BAM files.")

    def call_variants_per_read(self):
        """Call variants on each per-read BAM file using Minipileup."""
        self.create_folder(self.per_read_vcf_folder)
        bam_folders = [os.path.join(self.per_read_bam_folder, barcode) for barcode in self.barcode_list]
        with Pool(processes=self.threads) as pool:
            pool.map(self.process_bam_folder, bam_folders)
        print("Variant calling on per-read BAM files completed.")

    def process_bam_folder(self, barcode_bam_folder):
        """Process all per-read BAM files in a barcode folder."""
        barcode = os.path.basename(barcode_bam_folder)
        barcode_vcf_folder = os.path.join(self.per_read_vcf_folder, barcode)
        self.create_folder(barcode_vcf_folder)
        bam_files = [f for f in os.listdir(barcode_bam_folder) if f.endswith('.bam')]

        for bam_file in bam_files:
            bam_path = os.path.join(barcode_bam_folder, bam_file)
            vcf_file = bam_file.replace('.bam', '.vcf')
            vcf_path = os.path.join(barcode_vcf_folder, vcf_file)

            # Run Minipileup for each per-read BAM file
            minipileup_cmd = f'{self.minipileup_path} -vc -f {self.ref_genome} {bam_path} > {vcf_path}'
            try:
                subprocess.run(minipileup_cmd, shell=True, check=True)
                print(f"Called variants for {bam_file}, output saved to {vcf_file}")
            except subprocess.CalledProcessError as e:
                print(f"Minipileup failed for {bam_file}: {e}")

    def normalize_vcfs(self):
        """Normalize the per-read VCF files using bcftools norm."""
        self.create_folder(self.normalized_vcf_folder)
        vcf_folders = [os.path.join(self.per_read_vcf_folder, barcode) for barcode in self.barcode_list]
        with Pool(processes=self.threads) as pool:
            pool.map(self.process_vcf_folder, vcf_folders)
        print("Normalization of per-read VCF files completed.")

    def process_vcf_folder(self, barcode_vcf_folder):
        """Normalize all VCF files in a barcode folder."""
        barcode = os.path.basename(barcode_vcf_folder)
        barcode_norm_vcf_folder = os.path.join(self.normalized_vcf_folder, barcode)
        self.create_folder(barcode_norm_vcf_folder)
        vcf_files = [f for f in os.listdir(barcode_vcf_folder) if f.endswith('.vcf')]

        for vcf_file in vcf_files:
            input_vcf_path = os.path.join(barcode_vcf_folder, vcf_file)
            intermediate_vcf = os.path.join(barcode_norm_vcf_folder, f"{os.path.splitext(vcf_file)[0]}_norm_step1.vcf")
            final_vcf = os.path.join(barcode_norm_vcf_folder, f"{os.path.splitext(vcf_file)[0]}_norm_step2.vcf")

            # Run bcftools norm commands
            try:
                subprocess.run(['bcftools', 'norm', '-m-', '-o', intermediate_vcf, input_vcf_path], check=True)
                subprocess.run(['bcftools', 'norm', '-f', self.ref_genome, '-cw', '-o', final_vcf, intermediate_vcf], check=True)
                print(f"Normalized VCF file: {final_vcf}")
            except subprocess.CalledProcessError as e:
                print(f"bcftools norm failed for {vcf_file} in barcode {barcode}: {e}")

    def extract_variants_from_vcf(self, vcf_path, variant_names):
        """Extract variants from a VCF file and return the variant names present."""
        variants_present = []
        try:
            with open(vcf_path, 'r') as vcf_file:
                for line in vcf_file:
                    if line.startswith('#'):
                        continue
                    fields = line.strip().split('\t')
                    chrom = fields[0]
                    if chrom != self.gene_name:
                        continue  # Only process the selected gene
                    pos = int(fields[1])
                    ref = fields[3]
                    alt = fields[4]
                    # No need to split alt, since VCF is normalized
                    key = (chrom, pos, ref, alt)
                    if key in variant_names:
                        variants_present.append(variant_names[key])
        except Exception as e:
            print(f"Error reading VCF file {vcf_path}: {e}")
        return variants_present

    def collect_variants_and_build_matrix(self):
        """Collect variants from normalized VCF files and build co-occurrence matrices per barcode."""
        variant_counts_per_barcode = self.parse_known_variants_file()
        variant_freqs_per_barcode = self.parse_known_variants_file_frequencies()
        variants_list = list(set(self.variant_names.values()))
        variant_indices = {v: i for i, v in enumerate(variants_list)}

        for barcode in self.barcode_list:
            barcode_norm_vcf_folder = os.path.join(self.normalized_vcf_folder, barcode)
            vcf_files = [f for f in os.listdir(barcode_norm_vcf_folder) if f.endswith('_norm_step2.vcf')]
            cooccurrence_matrix = np.zeros((len(variants_list), len(variants_list)), dtype=int)
            total_reads_processed = 0

            # Set diagonal entries from known_variants_file (counts)
            variant_counts = variant_counts_per_barcode.get(barcode, {})
            for variant_name, count in variant_counts.items():
                if variant_name in variant_indices:
                    idx = variant_indices[variant_name]
                    cooccurrence_matrix[idx, idx] = count

            # Process per-read VCF files to fill off-diagonal entries
            for vcf_file in vcf_files:
                vcf_path = os.path.join(barcode_norm_vcf_folder, vcf_file)
                read_variants = self.extract_variants_from_vcf(vcf_path, self.variant_names)
                total_reads_processed += 1

                if len(read_variants) > 1:
                    # Update co-occurrence counts
                    for i in range(len(read_variants)):
                        vi = variant_indices[read_variants[i]]
                        for j in range(i + 1, len(read_variants)):
                            vj = variant_indices[read_variants[j]]
                            cooccurrence_matrix[vi, vj] += 1
                            cooccurrence_matrix[vj, vi] += 1  # Symmetric

            # Store total reads per barcode
            self.total_reads_per_barcode[barcode] = total_reads_processed

            # Save counts matrix
            cooccurrence_df = pd.DataFrame(cooccurrence_matrix, index=variants_list, columns=variants_list)
            counts_matrix_path = os.path.join(self.counts_matrices_folder, f"{barcode}_cooccurrence_matrix_counts.csv")
            cooccurrence_df.to_csv(counts_matrix_path)

            # Build frequency matrix
            frequencies_matrix = np.zeros_like(cooccurrence_matrix, dtype=float)

            # Total reads per barcode
            total_reads_per_barcode = total_reads_processed

            # Set diagonal entries from known_variants_file (frequencies)
            variant_freqs = variant_freqs_per_barcode.get(barcode, {})
            for variant_name, freq in variant_freqs.items():
                if variant_name in variant_indices:
                    idx = variant_indices[variant_name]
                    frequencies_matrix[idx, idx] = freq

            # Calculate off-diagonal frequencies
            for i in range(len(variants_list)):
                for j in range(len(variants_list)):
                    if i != j and total_reads_per_barcode > 0:
                        frequencies_matrix[i, j] = (cooccurrence_matrix[i, j] / total_reads_per_barcode) * 100

            # Round frequencies to two decimal places
            frequencies_matrix = np.round(frequencies_matrix, 2)

            # Save frequencies matrix
            frequencies_df = pd.DataFrame(frequencies_matrix, index=variants_list, columns=variants_list)
            frequencies_matrix_path = os.path.join(self.frequencies_matrices_folder, f"{barcode}_cooccurrence_matrix_frequencies.csv")
            frequencies_df.to_csv(frequencies_matrix_path)

            print(f"Counts and frequencies matrices saved for barcode {barcode}")

    def calculate_wildtype_frequencies_from_vcfs(self):
        """Calculate wild-type frequencies directly from VCF files based on total_reads_processed and save results to a file."""
        output_file = os.path.join(self.counts_matrices_folder, "wildtype_frequencies.csv")
        results = []

        for barcode in self.barcode_list:
            barcode_norm_vcf_folder = os.path.join(self.normalized_vcf_folder, barcode)
            if not os.path.exists(barcode_norm_vcf_folder):
                print(f"Normalized VCF folder not found for barcode {barcode}, skipping.")
                continue

            vcf_files = [f for f in os.listdir(barcode_norm_vcf_folder) if f.endswith('_norm_step2.vcf')]
            wildtype_count = 0

            for vcf_file in vcf_files:
                vcf_path = os.path.join(barcode_norm_vcf_folder, vcf_file)
                try:
                    with open(vcf_path, 'r') as vcf:
                        has_variant = False
                        for line in vcf:
                            if line.startswith('#'):
                                continue
                            fields = line.strip().split('\t')
                            chrom = fields[0]
                            if chrom != self.gene_name:
                                continue  # Only process the selected gene
                            pos = int(fields[1])
                            ref = fields[3]
                            alt = fields[4]
                            key = (chrom, pos, ref, alt)
                            if key in self.variant_names:
                                has_variant = True
                                break  # Variant found, not wildtype
                        if not has_variant:
                            wildtype_count += 1
                except Exception as e:
                    print(f"Error reading VCF file {vcf_path}: {e}")

            # Retrieve the total number of reads processed for this barcode
            total_reads = self.total_reads_per_barcode.get(barcode, 0)

            # Calculate wildtype frequency
            frequency_wildtype = round((wildtype_count / total_reads) * 100, 2) if total_reads > 0 else 0.0

            # Append results for this barcode
            results.append([barcode, wildtype_count, frequency_wildtype])

            # Debugging output
            print(f"Barcode: {barcode}")
            print(f"Total reads processed: {total_reads}")
            print(f"Wildtype reads: {wildtype_count}")
            print(f"Wildtype Frequency: {frequency_wildtype}%\n")

        # Create a DataFrame for the results
        results_df = pd.DataFrame(results, columns=['Sample', 'Count', 'Frequency'])

        # Function to extract the numeric part of the barcode for sorting
        def extract_barcode_number(sample_name):
            match = re.search(r'barcode(\d+)', sample_name)
            if match:
                return int(match.group(1))  # Extract the numeric part of the barcode
            return float('inf')  # If no barcode number found, place it at the end

        # Sort the DataFrame by barcode number
        results_df['Barcode_Num'] = results_df['Sample'].apply(extract_barcode_number)
        results_df = results_df.sort_values('Barcode_Num').drop('Barcode_Num', axis=1)

        # Write the results to a CSV file
        results_df.to_csv(output_file, index=False)

        print(f"Wild-type frequencies saved to {output_file}")

    def clean_up(self):
        """Delete intermediate directories after final output is generated."""
        dirs_to_remove = [self.per_read_bam_folder, self.per_read_vcf_folder, self.normalized_vcf_folder]
        for dir_path in dirs_to_remove:
            if os.path.exists(dir_path):
                shutil.rmtree(dir_path)
                print(f"Deleted directory: {dir_path}")
            else:
                print(f"Directory not found, skipping deletion: {dir_path}")

# Define gene positions
gene_positions = {
    'CYP51': [
        {'CHROM': 'CYP51', 'POS': 149, 'REF': 'T', 'ALT': 'C'},
        {'CHROM': 'CYP51', 'POS': 320, 'REF': 'A', 'ALT': 'T'},
        {'CHROM': 'CYP51', 'POS': 401, 'REF': 'A', 'ALT': 'G'},
        {'CHROM': 'CYP51', 'POS': 407, 'REF': 'T', 'ALT': 'C'},
        {'CHROM': 'CYP51', 'POS': 410, 'REF': 'A', 'ALT': 'T'},
        {'CHROM': 'CYP51', 'POS': 433, 'REF': 'A', 'ALT': 'C'},
        {'CHROM': 'CYP51', 'POS': 433, 'REF': 'A', 'ALT': 'T'},
        {'CHROM': 'CYP51', 'POS': 590, 'REF': 'A', 'ALT': 'G'},
        {'CHROM': 'CYP51', 'POS': 620, 'REF': 'G', 'ALT': 'A'},
        {'CHROM': 'CYP51', 'POS': 679, 'REF': 'T', 'ALT': 'A'},
        {'CHROM': 'CYP51', 'POS': 841, 'REF': 'C', 'ALT': 'T'},
        {'CHROM': 'CYP51', 'POS': 1398, 'REF': 'T', 'ALT': 'C'},
        {'CHROM': 'CYP51', 'POS': 1408, 'REF': 'C', 'ALT': 'G'},
        {'CHROM': 'CYP51', 'POS': 1413, 'REF': 'A', 'ALT': 'G'},
        {'CHROM': 'CYP51', 'POS': 1646, 'REF': 'CTAT', 'ALT': 'C'},
        {'CHROM': 'CYP51', 'POS': 1648, 'REF': 'A', 'ALT': 'C'},
        {'CHROM': 'CYP51', 'POS': 1649, 'REF': 'TGGC', 'ALT': 'T'},
        {'CHROM': 'CYP51', 'POS': 1653, 'REF': 'T', 'ALT': 'C'},
        {'CHROM': 'CYP51', 'POS': 1654, 'REF': 'A', 'ALT': 'C'},
        {'CHROM': 'CYP51', 'POS': 1811, 'REF': 'C', 'ALT': 'A'},
        {'CHROM': 'CYP51', 'POS': 1811, 'REF': 'C', 'ALT': 'G'},
        {'CHROM': 'CYP51', 'POS': 1854, 'REF': 'T', 'ALT': 'A'}
    ],
    'CytB': [
        {'CHROM': 'CytB', 'POS': 110, 'REF': 'G', 'ALT': 'T'},
        {'CHROM': 'CytB', 'POS': 428, 'REF': 'G', 'ALT': 'C'}
    ],
    'SdhB': [
        {'CHROM': 'SdhB', 'POS': 788, 'REF': 'A', 'ALT': 'C'},
        {'CHROM': 'SdhB', 'POS': 788, 'REF': 'A', 'ALT': 'T'},
        {'CHROM': 'SdhB', 'POS': 908, 'REF': 'G', 'ALT': 'C'},
        {'CHROM': 'SdhB', 'POS': 914, 'REF': 'A', 'ALT': 'T'},
        {'CHROM': 'SdhB', 'POS': 917, 'REF': 'C', 'ALT': 'T'},
        {'CHROM': 'SdhB', 'POS': 916, 'REF': 'A', 'ALT': 'G'}
    ],
    'SdhC': [
        {'CHROM': 'SdhC', 'POS': 365, 'REF': 'C', 'ALT': 'A'},
        {'CHROM': 'SdhC', 'POS': 368, 'REF': 'G', 'ALT': 'C'},
        {'CHROM': 'SdhC', 'POS': 379, 'REF': 'G', 'ALT': 'T'},
        {'CHROM': 'SdhC', 'POS': 386, 'REF': 'A', 'ALT': 'G'},
        {'CHROM': 'SdhC', 'POS': 387, 'REF': 'C', 'ALT': 'A'},
        {'CHROM': 'SdhC', 'POS': 387, 'REF': 'C', 'ALT': 'G'},
        {'CHROM': 'SdhC', 'POS': 508, 'REF': 'C', 'ALT': 'G'},
        {'CHROM': 'SdhC', 'POS': 582, 'REF': 'G', 'ALT': 'C'},
        {'CHROM': 'SdhC', 'POS': 582, 'REF': 'G', 'ALT': 'T'},
        {'CHROM': 'SdhC', 'POS': 581, 'REF': 'G', 'ALT': 'C'},
        {'CHROM': 'SdhC', 'POS': 581, 'REF': 'G', 'ALT': 'T'},
        {'CHROM': 'SdhC', 'POS': 580, 'REF': 'A', 'ALT': 'G'},
        {'CHROM': 'SdhC', 'POS': 584, 'REF': 'A', 'ALT': 'G'},
        {'CHROM': 'SdhC', 'POS': 611, 'REF': 'T', 'ALT': 'G'},
        {'CHROM': 'SdhC', 'POS': 625, 'REF': 'G', 'ALT': 'A'},
        {'CHROM': 'SdhC', 'POS': 632, 'REF': 'C', 'ALT': 'G'}
    ],
    'SdhD': [
        {'CHROM': 'SdhD', 'POS': 139, 'REF': 'C', 'ALT': 'T'},
        {'CHROM': 'SdhD', 'POS': 148, 'REF': 'A', 'ALT': 'T'},
        {'CHROM': 'SdhD', 'POS': 148, 'REF': 'A', 'ALT': 'C'},
        {'CHROM': 'SdhD', 'POS': 399, 'REF': 'A', 'ALT': 'G'},
        {'CHROM': 'SdhD', 'POS': 497, 'REF': 'A', 'ALT': 'G'},
        {'CHROM': 'SdhD', 'POS': 498, 'REF': 'C', 'ALT': 'A'},
        {'CHROM': 'SdhD', 'POS': 498, 'REF': 'C', 'ALT': 'G'}
    ]
}

# Set your paths and parameters here
bam_folder = sys.argv[1]
ref_genome = sys.argv[2]
gene_name = sys.argv[3]
minipileup_path = sys.argv[4]
output_folder = sys.argv[5]
counts_matrices_folder = sys.argv[6]
frequencies_matrices_folder = sys.argv[7]
known_variants_file = sys.argv[8]
threads = 30  # Set the number of threads for parallel processing

# Instantiate and run the workflow
workflow = VCFWorkflow(
    bam_folder=bam_folder,
    ref_genome=ref_genome,
    gene_name=gene_name,
    gene_positions=gene_positions.get(gene_name, []),
    minipileup_path=minipileup_path,
    output_folder=output_folder,
    counts_matrices_folder=counts_matrices_folder,
    frequencies_matrices_folder=frequencies_matrices_folder,
    known_variants_file=known_variants_file,
    threads=threads
)

# Run the workflow methods
workflow.split_bam_by_read()
workflow.call_variants_per_read()
workflow.normalize_vcfs()
workflow.collect_variants_and_build_matrix()
workflow.calculate_wildtype_frequencies_from_vcfs()

# Delete intermediate directories after final output is generated
workflow.clean_up()
