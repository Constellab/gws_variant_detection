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
                 counts_matrices_folder, frequencies_matrices_folder, known_variants_file, reads_to_process, threads):
        self.bam_folder = bam_folder
        self.ref_genome = ref_genome
        self.gene_name = gene_name
        self.gene_positions = gene_positions  # List of positions of interest
        self.minipileup_path = minipileup_path
        self.output_folder = output_folder
        self.counts_matrices_folder = counts_matrices_folder
        self.frequencies_matrices_folder = frequencies_matrices_folder
        self.known_variants_file = known_variants_file
        self.reads_to_process = reads_to_process  # User-provided total reads for wild-type calculation
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
            ('Cytb', 919, 'C', 'T'): 'S34L',
            ('Cytb', 1205, 'T', 'A'): 'F129L_A',
            ('Cytb', 1205, 'T', 'G'): 'F129L_G',
            ('Cytb', 1246, 'G', 'C'): 'G143A',
            ('Cytb', 1420, 'T', 'C'): 'L201S',
            ('Cytb', 1422, 'C', 'CATGAGG'): 'E203-DE-V204',
            ('Cytb', 1424, 'T', 'TGAGGTA'): 'E203-VE-V204',
            ('OSBP', 4169, 'T', 'G'): 'L863W',
            ('OSBP', 4091, 'A', 'T'): 'N837I',
            ('OSBP', 3890, 'G', 'T'): 'G770V',
            ('Tubulin', 1541, 'T', 'A'): 'C239S_AGC',
            ('Tubulin', 1541, 'T', 'G'): 'C239G_GGC',
            ('Tubulin', 1542, 'G', 'C'): 'C239S_TCC',
            ('ATPase', 2857, 'G', 'A'): 'G759S',
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
            subprocess.run(minipileup_cmd, shell=True)

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
            subprocess.run(['bcftools', 'norm', '-m-', '-o', intermediate_vcf, input_vcf_path], check=True)
            subprocess.run(['bcftools', 'norm', '-f', self.ref_genome, '-cw', '-o', final_vcf, intermediate_vcf], check=True)

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

    def calculate_wildtype_frequencies_from_matrices(self):
        """Calculate wild-type frequencies and save results to a file."""
        output_file = os.path.join(self.counts_matrices_folder, "wildtype_frequencies.csv")
        results = []

        for matrix_file in os.listdir(self.counts_matrices_folder):
            if matrix_file.endswith('_cooccurrence_matrix_counts.csv'):
                matrix_path = os.path.join(self.counts_matrices_folder, matrix_file)

                # Load co-occurrence matrix into DataFrame
                df_co_occ = pd.read_csv(matrix_path, index_col=0)

                # Convert the DataFrame to a numpy array for computation
                co_occurrence_matrix = np.array(df_co_occ)

                # Extract sample name from the file name
                sample_name = matrix_file.replace('_cooccurrence_matrix_counts.csv', '')

                # Total number of reads (user-provided)
                reads_to_process = self.reads_to_process

                # Sum of diagonal elements: reads with at least 1 variant
                reads_with_at_least_one_variant = np.sum(np.diag(co_occurrence_matrix))

                # Sum of the upper triangular part (excluding the diagonal): reads with exactly 2 variants
                upper_triangular_without_diagonal = np.triu(co_occurrence_matrix, k=1)
                reads_with_exactly_two_variants = np.sum(upper_triangular_without_diagonal)

                # Subtract the off-diagonal (reads with exactly 2 variants) to avoid double-counting
                reads_with_variants = reads_with_at_least_one_variant - reads_with_exactly_two_variants

                # Reads with 0 variants
                reads_with_no_variants = reads_to_process - reads_with_variants

                # If reads_with_no_variants is negative, set it to 0
                if reads_with_no_variants < 0:
                    reads_with_no_variants = 0

                # Calculate frequency using user-provided reads_to_process
                frequency_wildtype = round((reads_with_no_variants / reads_to_process) * 100, 2)

                # If frequency_wildtype is negative, set it to 0
                if frequency_wildtype < 0:
                    frequency_wildtype = 0

                # Append results for this sample
                results.append([sample_name, reads_with_no_variants, frequency_wildtype])

        # Create a DataFrame for results with updated column names
        results_df = pd.DataFrame(results, columns=['Sample', 'Count', 'Frequency'])

        # Function to extract the numeric barcode from the sample name
        def extract_barcode_number(sample_name):
            match = re.search(r'barcode(\d+)', sample_name)
            if match:
                return int(match.group(1))  # Extract the numeric barcode part as an integer
            return float('inf')  # If no barcode is found, put it at the end when sorting

        # Sort the DataFrame by the extracted barcode number
        results_df['Barcode_Num'] = results_df['Sample'].apply(extract_barcode_number)
        results_df = results_df.sort_values('Barcode_Num').drop('Barcode_Num', axis=1)

        # Write sorted DataFrame to a CSV file
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

# Set your paths and parameters here
bam_folder = sys.argv[1]
ref_genome = sys.argv[2]
gene_name = sys.argv[3]
minipileup_path = sys.argv[4]
output_folder = sys.argv[5]
counts_matrices_folder = sys.argv[6]
frequencies_matrices_folder = sys.argv[7]
known_variants_file = sys.argv[8]
reads_to_process = int(sys.argv[9])
threads = 30  # Set the number of threads for parallel processing

# Define gene positions
gene_positions = {
    'Cytb': [
        {'CHROM': 'Cytb', 'POS': 919, 'REF': 'C', 'ALT': 'T'},
        {'CHROM': 'Cytb', 'POS': 1205, 'REF': 'T', 'ALT': 'A'},
        {'CHROM': 'Cytb', 'POS': 1205, 'REF': 'T', 'ALT': 'G'},
        {'CHROM': 'Cytb', 'POS': 1246, 'REF': 'G', 'ALT': 'C'},
        {'CHROM': 'Cytb', 'POS': 1420, 'REF': 'T', 'ALT': 'C'},
        {'CHROM': 'Cytb', 'POS': 1422, 'REF': 'C', 'ALT': 'CATGAGG'},
        {'CHROM': 'Cytb', 'POS': 1424, 'REF': 'T', 'ALT': 'TGAGGTA'}
    ],
    'OSBP': [
        {'CHROM': 'OSBP', 'POS': 4169, 'REF': 'T', 'ALT': 'G'},
        {'CHROM': 'OSBP', 'POS': 4091, 'REF': 'A', 'ALT': 'T'},
        {'CHROM': 'OSBP', 'POS': 3890, 'REF': 'G', 'ALT': 'T'}
     ],
    'Tubulin': [
        {'CHROM': 'Tubulin', 'POS': 1541, 'REF': 'T', 'ALT': 'A'},
        {'CHROM': 'Tubulin', 'POS': 1541, 'REF': 'T', 'ALT': 'G'},
        {'CHROM': 'Tubulin', 'POS': 1542, 'REF': 'G', 'ALT': 'C'}
    ],
    'ATPase': [
        {'CHROM': 'ATPase', 'POS': 2857, 'REF': 'G', 'ALT': 'A'}
    ]
}

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
    reads_to_process=reads_to_process,
    threads=threads
)

# Run the workflow methods
workflow.split_bam_by_read()
workflow.call_variants_per_read()
workflow.normalize_vcfs()
workflow.collect_variants_and_build_matrix()
workflow.calculate_wildtype_frequencies_from_matrices()

# Delete intermediate directories after final output is generated
workflow.clean_up()
