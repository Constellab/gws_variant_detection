import os
import pandas as pd
import subprocess
import numpy as np
import re
import sys

class VCFWorkflow:
    def __init__(self, input_vcf_folder, input_bam_folder, reference_genome, gene_name, intermediate_folder, output_folder, output_matrix_folder, freq_matrix_folder, known_variants_file, reads_to_process):
        self.input_vcf_folder = input_vcf_folder
        self.input_bam_folder = input_bam_folder
        self.reference_genome = reference_genome
        self.gene_name = gene_name
        self.intermediate_folder = intermediate_folder
        self.output_folder = output_folder
        self.output_matrix_folder = output_matrix_folder  # Folder for count-based matrices
        self.freq_matrix_folder = freq_matrix_folder  # Folder for frequency-based matrices
        self.known_variants_file = known_variants_file
        self.gene_positions = self.get_gene_positions()
        self.reads_to_process = reads_to_process  # Pass total reads 5000

    def create_folder(self, folder_path):
        """Ensure a folder exists; if not, create it."""
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            print(f"Created folder: {folder_path}")

    def get_gene_positions(self):
        """Retrieve predefined positions of interest for the specified gene."""
        positions = {
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
        return positions.get(self.gene_name, [])

    def extract_variants(self):
        """Extracts variants based on gene-specific positions."""
        self.create_folder(self.intermediate_folder)
        for vcf_file in os.listdir(self.input_vcf_folder):
            if vcf_file.endswith('.vcf') or vcf_file.endswith('.vcf.gz'):
                input_path = os.path.join(self.input_vcf_folder, vcf_file)
                output_path = os.path.join(self.intermediate_folder, vcf_file.replace('.gz', ''))

                # Read VCF into DataFrame, ensuring header is maintained
                df = pd.read_csv(input_path, comment='#', sep='\t', dtype={'POS': int},
                                 names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'],
                                 compression='gzip' if input_path.endswith('.gz') else None)

                # Filter DataFrame for required positions
                filtered_df = df[df.apply(lambda row: any((row['CHROM'] == pos['CHROM'] and row['POS'] == pos['POS'] and row['REF'] == pos['REF'] and row['ALT'] == pos['ALT']) for pos in self.gene_positions), axis=1)]

                # Write filtered DataFrame to VCF, preserving the original headers
                with open(output_path, 'w') as f:
                    # Re-include the VCF headers, typically found in the original file
                    with open(input_path, 'r') as original_vcf:
                        for line in original_vcf:
                            if line.startswith('#'):
                                f.write(line)
                    filtered_df.to_csv(f, sep='\t', index=False, header=False)

    def phase_vcfs(self):
        """Phase the extracted VCF files using corresponding BAM files."""
        self.create_folder(self.output_folder)
        self.create_folder(self.output_matrix_folder)  # Create the output matrix folder
        self.create_folder(self.freq_matrix_folder)  # Create the frequency matrix folder
        for vcf_file in os.listdir(self.intermediate_folder):
            if vcf_file.endswith('.vcf'):
                vcf_prefix = vcf_file.split('-')[0]
                bam_file = self.find_bam_file(vcf_prefix)
                if bam_file:
                    input_vcf = os.path.join(self.intermediate_folder, vcf_file)
                    input_bam = os.path.join(self.input_bam_folder, bam_file)
                    output_vcf = os.path.join(self.output_folder, f"{vcf_file}")
                    self.phase_vcf(input_vcf, input_bam, output_vcf)
                    self.assign_variant_names(output_vcf, vcf_prefix)  # Call the variant name assignment function
                else:
                    print(f"No corresponding BAM file found for {vcf_file}")

    def find_bam_file(self, vcf_prefix):
        """Find the corresponding BAM file for a given VCF file."""
        for bam_file in os.listdir(self.input_bam_folder):
            if bam_file.startswith(vcf_prefix) and bam_file.endswith('.bam'):
                return bam_file
        return None

    def phase_vcf(self, input_vcf, input_bam, output_vcf):
        """Run the phasing command using Whatshap."""
        command = [
            'whatshap', 'phase', '--ignore-read-groups', '--mapping-quality', '0', '--default-gq', '0',
            '-o', output_vcf,
            '--reference', self.reference_genome,
            input_vcf, input_bam
        ]
        subprocess.run(command, check=True)

    def assign_variant_names(self, output_vcf, barcode_prefix):
        """Assign variant names based on CHROM, POS, REF, ALT."""
        # Read the phased VCF into DataFrame
        df = pd.read_csv(output_vcf, comment='#', sep='\t', dtype={'POS': int},
                         names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'])

        # Mapping of variant names based on CHROM, POS, REF, ALT
        variant_names = {
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

        # Assign variant names based on CHROM, POS, REF, ALT
        df['variant_name'] = df.apply(lambda row: variant_names.get((row['CHROM'], row['POS'], row['REF'], row['ALT']), 'Unknown'), axis=1)

        # Rearrange the columns so that variant_name is immediately after POS
        df = df[['CHROM', 'POS', 'variant_name', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']]

        # Write the DataFrame back to the VCF file
        df.to_csv(output_vcf, sep='\t', index=False)

        # Print the phased VCF file
        print("Phased VCF File:")
        print(df)

        # Indicate which variants should be phased together
        print("Phased Variants:")
        if 'PS' in df.columns:
            phased_variants = df.dropna(subset=['PS']).groupby('PS')
            for ps, group in phased_variants:
                variants = group['variant_name'].values
                print(f"Phase set {ps}: {variants}")
        else:
            print("No 'PS' column found in the VCF file.")

        # Parse the modified VCF and build co-occurrence matrices
        vcf_df = parse_vcf(output_vcf)
        cooccurrence_matrix_counts = self.build_cooccurrence_matrix(vcf_df, barcode_prefix, 'COUNT')
        cooccurrence_matrix_freqs = self.build_cooccurrence_matrix(vcf_df, barcode_prefix, 'FREQ')

        # Save the matrices to CSV files
        output_matrix_path_counts = os.path.join(self.output_matrix_folder, f"{os.path.basename(output_vcf).split('.')[0]}_matrix_counts.csv")
        output_matrix_path_freqs = os.path.join(self.freq_matrix_folder, f"{os.path.basename(output_vcf).split('.')[0]}_matrix_freqs.csv")
        cooccurrence_matrix_counts.to_csv(output_matrix_path_counts)
        cooccurrence_matrix_freqs.to_csv(output_matrix_path_freqs)
        print(f"Co-occurrence matrix (counts) saved to {output_matrix_path_counts}.")
        print(f"Co-occurrence matrix (freqs) saved to {output_matrix_path_freqs}.")

    def build_cooccurrence_matrix(self, vcf_df, barcode_prefix, value_type):
        # Read the known variants file
        known_variants_df = pd.read_csv(self.known_variants_file, sep='\t')

        # Extract variant names and their values (counts or frequencies) for the specific barcode
        variant_values = {}
        value_suffix = 'VARIANT_COUNT' if value_type == 'COUNT' else 'VARIANT_FREQ'
        for variant_name in known_variants_df['variant_name']:
            value_column = f"{barcode_prefix}_{value_suffix}"
            if value_column in known_variants_df.columns:
                value = known_variants_df.loc[known_variants_df['variant_name'] == variant_name, value_column].values[0]
                if value_type == 'FREQ':
                    value = float(value.strip('%'))  # Remove '%' and convert to float
                variant_values[variant_name] = value

        # Debugging output for variant values
        print("Variant values extracted from known variants file:")
        for variant, value in variant_values.items():
            print(f"{variant}: {value}")

        # All variants, including unphased ones
        all_variants = vcf_df['variant_name'].unique()
        cooccurrence_df = pd.DataFrame(0, index=all_variants, columns=all_variants, dtype=float if value_type == 'FREQ' else int)

        # Initially populate the diagonal with values for each variant using known_variants_df
        for variant in all_variants:
            if variant in variant_values:
                cooccurrence_df.at[variant, variant] = variant_values[variant]

        # Process phased variants for off-diagonal entries
        if 'PS' in vcf_df.columns:
            phased_variants = vcf_df.dropna(subset=['PS'])
            phase_groups = phased_variants.groupby('PS')

            for ps, group in phase_groups:
                variants = group['variant_name'].values
                print(f"Phase set {ps}: {variants}")  # Debugging output
                for i, variant_i in enumerate(variants):
                    for j, variant_j in enumerate(variants):
                        if i != j:
                            # Use the minimum value from known_variants_df for the phased variants
                            value_i = variant_values.get(variant_i, 0)
                            value_j = variant_values.get(variant_j, 0)
                            min_value = min(value_i, value_j)
                            print(f"Updating co-occurrence for {variant_i} and {variant_j} with value {min_value}")
                            cooccurrence_df.at[variant_i, variant_j] = min_value
                            cooccurrence_df.at[variant_j, variant_i] = min_value  # Ensure symmetry

                            # Debugging output to check if values are being updated
                            print(f"Updated co-occurrence for {variant_i} and {variant_j} with value {min_value}")

        # Debugging output to print the co-occurrence matrix
        print("Co-occurrence matrix:")
        print(cooccurrence_df)

        return cooccurrence_df

    def calculate_wildtype_frequencies_from_matrices(self):
        """Calculate wild-type frequencies and save results to a file."""
        output_file = os.path.join(self.output_matrix_folder, "wildtype_frequencies.csv")
        results = []

        for matrix_file in os.listdir(self.output_matrix_folder):
            if matrix_file.endswith('_matrix_counts.csv'):
                matrix_path = os.path.join(self.output_matrix_folder, matrix_file)

                # Load co-occurrence matrix into DataFrame
                df_co_occ = pd.read_csv(matrix_path, index_col=0)

                # Convert the DataFrame to a numpy array for computation
                co_occurrence_matrix = np.array(df_co_occ)

                # Total number of reads (defined when initializing the workflow)
                total_reads = self.reads_to_process

                # Sum of diagonal elements: reads with at least 1 variant
                reads_with_at_least_one_variant = np.sum(np.diag(co_occurrence_matrix))

                # Sum of the upper triangular part (excluding the diagonal): reads with exactly 2 variants
                upper_triangular_without_diagonal = np.triu(co_occurrence_matrix, k=1)
                reads_with_exactly_two_variants = np.sum(upper_triangular_without_diagonal)

                # Subtract the off-diagonal (reads with exactly 2 variants) to avoid double-counting
                reads_with_variants = reads_with_at_least_one_variant - reads_with_exactly_two_variants

                # Reads with 0 variants
                reads_with_no_variants = total_reads - reads_with_variants

                # If sum_no_mut (reads_with_no_variants) is negative, set it to 0
                if reads_with_no_variants < 0:
                    reads_with_no_variants = 0

                # Store only the "Reads with 0 variants" count and percentage
                frequency_wildtype = round((reads_with_no_variants / total_reads) * 100, 2)

                # If frequency_wildtype is negative, set it to 0
                if frequency_wildtype < 0:
                    frequency_wildtype = 0

                # Extract sample name from the file name
                sample_name = matrix_file.replace('_matrix_counts.csv', '')

                # Append results for this sample
                results.append([sample_name, reads_with_no_variants, frequency_wildtype])

        # Create a DataFrame for results
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

def parse_vcf(vcf_path):
    # Read VCF into DataFrame, skipping headers that start with '##'
    vcf_df = pd.read_csv(vcf_path, comment='#', sep='\t', header=None,
                         names=['CHROM', 'POS', 'variant_name', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'DATA'])

    # Remove the first line if it contains column names
    if vcf_df.iloc[0]['CHROM'] == 'CHROM' and vcf_df.iloc[0]['POS'] == 'POS':
        vcf_df = vcf_df.iloc[1:]

    # Split FORMAT and DATA fields into separate columns
    format_fields = vcf_df['FORMAT'].str.split(':', expand=True)
    data_fields = vcf_df['DATA'].apply(lambda x: pd.Series(x.split(':')))

    # Match FORMAT fields with data columns
    for i, field in enumerate(format_fields.iloc[0]):
        vcf_df[field] = data_fields[i].apply(lambda x: x.split(',') if isinstance(x, str) and ',' in x else x)

    # Extract AD and PS from the data, handling missing PS and AD fields
    if 'PS' in format_fields.iloc[0].values:
        vcf_df['PS'] = data_fields[format_fields.iloc[0].tolist().index('PS')].apply(pd.to_numeric, errors='coerce')
    else:
        vcf_df['PS'] = None

    vcf_df['AD'] = vcf_df['AD'].apply(lambda x: [int(y) for y in x] if isinstance(x, list) and x[0] != '.' else [0, 0])

    return vcf_df


input_vcf_folder = sys.argv[1]
input_bam_folder = sys.argv[2]
reference_genome = sys.argv[3]
gene_name = str(sys.argv[4])
intermediate_folder = sys.argv[5]
output_folder = sys.argv[6]
output_matrix_folder = sys.argv[7]
freq_matrix_folder = sys.argv[8]
known_variants_file = sys.argv[9]
reads_to_process = int(sys.argv[10])

workflow = VCFWorkflow(input_vcf_folder, input_bam_folder, reference_genome, gene_name, intermediate_folder, output_folder, output_matrix_folder, freq_matrix_folder, known_variants_file, reads_to_process)
workflow.extract_variants()
workflow.phase_vcfs()
workflow.calculate_wildtype_frequencies_from_matrices()


