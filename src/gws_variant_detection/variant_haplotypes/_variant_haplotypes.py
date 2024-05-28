import os
import pandas as pd
import subprocess
import sys

class VCFWorkflow:
    def __init__(self, input_vcf_folder, input_bam_folder, reference_genome, gene_name, intermediate_folder, output_folder, output_matrix_folder):
        self.input_vcf_folder = input_vcf_folder
        self.intermediate_folder = intermediate_folder
        self.input_bam_folder = input_bam_folder
        self.output_folder = output_folder
        self.output_matrix_folder = output_matrix_folder  # New folder for output matrices
        self.reference_genome = reference_genome
        self.gene_name = gene_name
        self.gene_positions = self.get_gene_positions()

    def create_folder(self, folder_path):
        """Ensure a folder exists; if not, create it."""
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            print(f"Created folder: {folder_path}")

    def get_gene_positions(self):
        """Retrieve predefined positions of interest for the specified gene."""
        positions = {
            'Cytb': [
                {'CHROM': 'Cytb', 'POS': 919, 'REF': 'C', 'ALT': 'T'},
                {'CHROM': 'Cytb', 'POS': 1205, 'REF': 'T', 'ALT': 'A'},
                {'CHROM': 'Cytb', 'POS': 1205, 'REF': 'T', 'ALT': 'G'},
                {'CHROM': 'Cytb', 'POS': 1246, 'REF': 'G', 'ALT': 'C'},
                {'CHROM': 'Cytb', 'POS': 1420, 'REF': 'T', 'ALT': 'C'},
                {'CHROM': 'Cytb', 'POS': 1422, 'REF': 'C', 'ALT': 'CATGAGG'},
                {'CHROM': 'Cytb', 'POS': 1431, 'REF': 'G', 'ALT': 'GAGGTAG'}
            ],
            'OSBP': [
                {'CHROM': 'OSBP', 'POS': 3889, 'REF': 'G', 'ALT': 'T'},
                {'CHROM': 'OSBP', 'POS': 4090, 'REF': 'A', 'ALT': 'T'},
                {'CHROM': 'OSBP', 'POS': 4090, 'REF': 'A', 'ALT': 'ATT'},
                {'CHROM': 'OSBP', 'POS': 4090, 'REF': 'A', 'ALT': 'ATC'},
                {'CHROM': 'OSBP', 'POS': 4168, 'REF': 'T', 'ALT': 'G'}
            ],
            'Tubulin': [
                {'CHROM': 'Tubulin', 'POS': 1541, 'REF': 'T', 'ALT': 'A'},
                {'CHROM': 'Tubulin', 'POS': 1541, 'REF': 'T', 'ALT': 'G'},
                {'CHROM': 'Tubulin', 'POS': 1542, 'REF': 'G', 'ALT': 'C'}
            ],
            'ATPase': []
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
        for vcf_file in os.listdir(self.intermediate_folder):
            if vcf_file.endswith('.vcf'):
                vcf_prefix = vcf_file.split('.')[0]
                bam_file = self.find_bam_file(vcf_prefix)
                if bam_file:
                    input_vcf = os.path.join(self.intermediate_folder, vcf_file)
                    input_bam = os.path.join(self.input_bam_folder, bam_file)
                    output_vcf = os.path.join(self.output_folder, f"{vcf_file}")
                    self.phase_vcf(input_vcf, input_bam, output_vcf)
                    self.assign_variant_names(output_vcf)  # Call the variant name assignment function
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

    def assign_variant_names(self, output_vcf):
        """Assign variant names based on CHROM, POS, REF, ALT."""
        # Read the phased VCF into DataFrame
        df = pd.read_csv(output_vcf, comment='#', sep='\t', dtype={'POS': int},
                         names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'])

        # Mapping of variant names based on CHROM, POS, REF, ALT
        variant_names = {
            ('Cytb', 919, 'C', 'T'): 'Variant4_S34L',
            ('Cytb', 1205, 'T', 'A'): 'Variant5_F129L_A',
            ('Cytb', 1205, 'T', 'G'): 'Variant5_F129L_G',
            ('Cytb', 1246, 'G', 'C'): 'Variant6_G143A',
            ('Cytb', 1420, 'T', 'C'): 'Variant3_L201S',
            ('Cytb', 1422, 'C', 'CATGAGG'): 'Variant1_E203-DE-V204',
            ('Cytb', 1431, 'G', 'GAGGTAG'): 'Variant2_E203-VE-V204',
            ('OSBP', 3889, 'G', 'T'): 'G770V',
            ('OSBP', 4090, 'A', 'T'): 'N837I',
            ('OSBP', 4090, 'A', 'ATT'): 'N837I_ATT',
            ('OSBP', 4090, 'A', 'ATC'): 'N837I_ATC',
            ('OSBP', 4168, 'T', 'G'): 'L863W',
            ('Tubulin', 1541, 'T', 'A'): 'C239S_AGC',
            ('Tubulin', 1541, 'T', 'G'): 'C239G_GGC',
            ('Tubulin', 1542, 'G', 'C'): 'C239S_TCC',
            # Add more mappings for other genes if needed
        }

        # Assign variant names based on CHROM, POS, REF, ALT
        df['variant_name'] = df.apply(lambda row: variant_names.get((row['CHROM'], row['POS'], row['REF'], row['ALT']), 'Unknown'), axis=1)

        # Rearrange the columns so that variant_name is immediately after POS
        df = df[['CHROM', 'POS', 'variant_name', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']]

        # Write the DataFrame back to the VCF file
        df.to_csv(output_vcf, sep='\t', index=False)

        # Parse the modified VCF and build co-occurrence matrix
        vcf_df = parse_vcf(output_vcf)
        cooccurrence_matrix = build_cooccurrence_matrix(vcf_df)

        # Save the matrix to a CSV file
        output_matrix_path = os.path.join(self.output_matrix_folder, f"{os.path.basename(output_vcf).split('.')[0]}_matrix.csv")
        cooccurrence_matrix.to_csv(output_matrix_path)
        print(f"Co-occurrence matrix saved to {output_matrix_path}.")

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
    if 'PS' in vcf_df.columns:
        vcf_df['PS'] = pd.to_numeric(vcf_df['PS'], errors='coerce')
    else:
        vcf_df['PS'] = None

    vcf_df['AD'] = vcf_df['AD'].apply(lambda x: [int(y) for y in x] if isinstance(x, list) and x[0] != '.' else [0, 0])

    return vcf_df

def build_cooccurrence_matrix(vcf_df):
    # All variants, including unphased ones
    all_variants = vcf_df['variant_name'].unique()
    cooccurrence_df = pd.DataFrame(0, index=all_variants, columns=all_variants, dtype=int)

    # Initially populate the diagonal with allele counts for each variant
    for _, row in vcf_df.iterrows():
        variant = row['variant_name']
        allele_count = row['AD'][1] if len(row['AD']) > 1 else 0
        cooccurrence_df.at[variant, variant] += allele_count

    # Process phased variants for off-diagonal entries
    phased_variants = vcf_df.dropna(subset=['PS'])
    phase_groups = phased_variants.groupby('PS')

    for _, group in phase_groups:
        variants = group['variant_name'].values
        ad_values = group['AD'].apply(lambda x: x[1] if len(x) > 1 else 0)

        for i, variant_i in enumerate(variants):
            for j, variant_j in enumerate(variants):
                if i != j:
                    min_value = min(ad_values.iloc[i], ad_values.iloc[j])
                    cooccurrence_df.at[variant_i, variant_j] += min_value

    return cooccurrence_df

# Example usage
input_vcf_folder = sys.argv[1]
input_bam_folder = sys.argv[2]
reference_genome = sys.argv[3]
gene_name = str(sys.argv[4])
intermediate_folder = sys.argv[5]
output_folder = sys.argv[6]
output_matrix_folder = sys.argv[7]

workflow = VCFWorkflow(input_vcf_folder, input_bam_folder, reference_genome, gene_name, intermediate_folder, output_folder, output_matrix_folder)
workflow.extract_variants()
workflow.phase_vcfs()
