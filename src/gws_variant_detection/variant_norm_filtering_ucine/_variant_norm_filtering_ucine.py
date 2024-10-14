import sys
import os
import subprocess
import pandas as pd

class VCFNormalizer:
    def __init__(
            self, vcf_folder, ref_genome, gene_name, reads_to_process, final_normalized_output_folder, output_folder,
            min_qual):
        self.vcf_folder = vcf_folder
        self.ref_genome = ref_genome
        self.output_folder = output_folder
        self.final_normalized_output_folder = final_normalized_output_folder
        self.gene_name = gene_name
        self.reads_to_process = reads_to_process
        self.min_qual = min_qual
        print("Type of reads_to_process:", type(self.reads_to_process))

        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        if not os.path.exists(self.final_normalized_output_folder):
            os.makedirs(self.final_normalized_output_folder)

    def normalize_vcf_files(self):
        for vcf_file in os.listdir(self.vcf_folder):
            if vcf_file.endswith('.vcf'):
                self.process_vcf_file(vcf_file)

    def process_vcf_file(self, vcf_file):
        final_vcf = os.path.join(self.final_normalized_output_folder, f"{os.path.splitext(vcf_file)[0]}_norm_step2.vcf")
        intermediate_vcf = os.path.join(self.output_folder, f"{os.path.splitext(vcf_file)[0]}_norm_step1.vcf")
        input_path = os.path.join(self.vcf_folder, vcf_file)
        subprocess.run(['bcftools', 'norm', '-m-', '-o', intermediate_vcf, input_path], check=True)
        subprocess.run(['bcftools', 'norm', '-f', self.ref_genome, '-cw',
                       '-o', final_vcf, intermediate_vcf], check=True)

        self.extract_data_to_tsv(final_vcf, vcf_file)

    def extract_data_to_tsv(self, final_vcf_path, vcf_file):
        barcode = os.path.splitext(vcf_file)[0]  # Assuming barcode is derived from the file name
        output_tsv = os.path.join(self.output_folder, f"{barcode}_variants.tsv")
        with open(final_vcf_path, 'r') as vcf, open(output_tsv, 'w') as tsv:
            tsv.write("CHROM\tPOS\tREF\tALT\t{}_AD\t{}_VARIANT_FREQ\t{}_VARIANT_COUNT\t{}_QUAL\n".format(
                barcode, barcode, barcode, barcode))
            for line in vcf:
                if not line.startswith('#'):
                    columns = line.strip().split('\t')
                    pos = int(columns[1])
                    qual = float(columns[5])
                    if qual >= self.min_qual:
                        format_parts = columns[8].split(':')
                        ad_index = format_parts.index('AD') if 'AD' in format_parts else None
                        if ad_index is not None and len(columns) > 9:
                            sample_info = columns[9].split(':')
                            ad_value = sample_info[ad_index]
                            ad_parts = ad_value.split(',')
                            if len(ad_parts) == 2 and all(x.isdigit() for x in ad_parts):
                                ref_count = int(ad_parts[0])
                                alt_count = int(ad_parts[1])
                                total_count = ref_count + alt_count
                                variant_freq = (alt_count / total_count * 100) if total_count > 0 else 0
                                variant_count = round((variant_freq / 100) * self.reads_to_process)
                                tsv.write(
                                    f"{columns[0]}\t{pos}\t{columns[3]}\t{columns[4]}\t{ad_value}\t{variant_freq:.2f}%\t{variant_count}\t{qual}\n")

    def merge_variant_files(self):
        tsv_files = [f for f in os.listdir(self.output_folder) if f.endswith('_variants.tsv')]
        merged_df = None
        for tsv in tsv_files:
            df = pd.read_csv(os.path.join(self.output_folder, tsv), sep='\t')
            barcode = tsv.split('_variants.tsv')[0]  # Extract barcode from file name

            # Rename columns to include barcode in the column name
            df.rename(columns={
                'AD': f'{barcode}_AD',
                'VARIANT_FREQ': f'{barcode}_VARIANT_FREQ',
                'VARIANT_COUNT': f'{barcode}_VARIANT_COUNT',
                'QUAL': f'{barcode}_QUAL'
            }, inplace=True)

            if merged_df is None:
                merged_df = df
            else:
                merged_df = pd.merge(merged_df, df, on=['CHROM', 'POS', 'REF', 'ALT'], how='outer')

        merged_df.fillna(0, inplace=True)  # Fill missing values with 0
        merged_df.to_csv(os.path.join(self.output_folder, 'merged_variants.tsv'), sep='\t', index=False)
        return merged_df

    def parse_variants(self, merged_df):
        # Define known variants for each gene including new ones
        gene_variants = {
            'CYP51': [
                {'variant_name': 'Y136F', 'CHROM': 'CYP51', 'POS': 407, 'REF': 'A', 'ALT': 'T'}
            ],
            'CytB': [
                {'variant_name': 'G143A', 'CHROM': 'CytB', 'POS': 416, 'REF': 'G', 'ALT': 'C'}
            ],
            'SdhB': [
                {'variant_name': 'H242Y', 'CHROM': 'SdhB', 'POS': 1007, 'REF': 'C', 'ALT': 'T'},
                {'variant_name': 'H242R', 'CHROM': 'SdhB', 'POS': 1008, 'REF': 'A', 'ALT': 'G'},
                {'variant_name': 'H242L', 'CHROM': 'SdhB', 'POS': 1008, 'REF': 'A', 'ALT': 'T'},
                {'variant_name': 'I244V', 'CHROM': 'SdhB', 'POS': 1013, 'REF': 'A', 'ALT': 'G'}
            ],
            'SdhC': [
                {'variant_name': 'G25R', 'CHROM': 'SdhC', 'POS': 249, 'REF': 'G', 'ALT': 'A'},
                {'variant_name': 'A83V', 'CHROM': 'SdhC', 'POS': 543, 'REF': 'C', 'ALT': 'T'},
                {'variant_name': 'T98I', 'CHROM': 'SdhC', 'POS': 588, 'REF': 'C', 'ALT': 'T'},
                {'variant_name': 'A128T', 'CHROM': 'SdhC', 'POS': 677, 'REF': 'G', 'ALT': 'A'},
                {'variant_name': 'D155N', 'CHROM': 'SdhC', 'POS': 758, 'REF': 'G', 'ALT': 'A'},
                {'variant_name': 'G169S', 'CHROM': 'SdhC', 'POS': 800, 'REF': 'G', 'ALT': 'A'},
                {'variant_name': 'G169D', 'CHROM': 'SdhC', 'POS': 801, 'REF': 'G', 'ALT': 'A'}
            ],
            'SdhD': [
                {'variant_name': 'K37Q', 'CHROM': 'SdhD', 'POS': 185, 'REF': 'A', 'ALT': 'C'},
                {'variant_name': 'I71F', 'CHROM': 'SdhD', 'POS': 287, 'REF': 'A', 'ALT': 'T'},
                {'variant_name': 'N69D', 'CHROM': 'SdhD', 'POS': 281, 'REF': 'A', 'ALT': 'G'},
                {'variant_name': 'I136F', 'CHROM': 'SdhD', 'POS': 588, 'REF': 'A', 'ALT': 'T'}
            ]
        }

        known_variants = gene_variants.get(self.gene_name, [])
        known_df = pd.DataFrame(known_variants)
        merged_known = pd.merge(known_df, merged_df, on=['CHROM', 'POS', 'REF', 'ALT'], how='left')

        for col in merged_known.columns:
            if '_VARIANT_FREQ' in col:
                merged_known[col] = merged_known[col].fillna('0%')
            elif '_AD' in col or '_VARIANT_COUNT' in col:
                merged_known[col] = merged_known[col].fillna(0)

        for col in merged_known.columns:
            if '_QUAL' in col:
                merged_known.drop(col, axis=1, inplace=True)

        # Save known variants
        merged_known.to_csv(os.path.join(self.output_folder, 'known_variants.tsv'), sep='\t', index=False)

        # Find unknown variants
        unknown = merged_df[~merged_df.set_index(['CHROM', 'POS', 'REF', 'ALT']).index.isin(
            known_df.set_index(['CHROM', 'POS', 'REF', 'ALT']).index)]

        for col in unknown.columns:
            if '_QUAL' in col:
                unknown.drop(col, axis=1, inplace=True)

        # Save unknown variants
        unknown.to_csv(os.path.join(self.output_folder, 'unknown_variants.tsv'), sep='\t', index=False)



vcf_folder = sys.argv[1]
ref_genome = sys.argv[2]
gene_name = str(sys.argv[3])  # Modify as needed for the gene of interest
reads_to_process = int(sys.argv[4])
final_normalized_output_folder = sys.argv[5]
output_folder = sys.argv[6]
min_qual = 0

normalizer = VCFNormalizer(vcf_folder, ref_genome, gene_name, reads_to_process,
                           final_normalized_output_folder, output_folder, min_qual)
normalizer.normalize_vcf_files()
merged_df = normalizer.merge_variant_files()
normalizer.parse_variants(merged_df)

# Cleanup
for file in os.listdir(output_folder):
    if file.endswith(".vcf"):
        os.remove(os.path.join(output_folder, file))

# Function to properly rename the columns


# Function to properly rename the columns
def rename_columns(df):
    new_columns = []
    for col in df.columns:
        if col == 'variant_name':  # Keep the 'variant_name' column unchanged
            new_columns.append(col)
        else:
            # Split by hyphen to extract the barcode part before '-'
            barcode = col.split('-')[0]
            # Split by underscore to extract the method suffix (like VARIANT_COUNT, VARIANT_FREQ)
            suffix = col.split('_')[-1]

            # Ensure that the suffix is correctly named as VARIANT_COUNT or VARIANT_FREQ
            if "COUNT" in col:
                suffix = "VARIANT_COUNT"
            elif "FREQ" in col:
                suffix = "VARIANT_FREQ"

            # Combine the barcode with the corrected suffix
            new_column_name = f"{barcode}_{suffix}"
            new_columns.append(new_column_name)

    # Assign the newly formed column names to the dataframe
    df.columns = new_columns
    return df

# Function to process the known_variants DataFrame
def modify_dataframe(df):
    # Rename the columns
    df = rename_columns(df)

    # Sort columns by barcode and method (VARIANT_COUNT/VARIANT_FREQ)
    sorted_columns = sorted(df.columns, key=lambda x: (x.split('_')[0], x.split('_')[-1]))

    # Reorder DataFrame columns with 'variant_name' first
    df = df[['variant_name'] + [col for col in sorted_columns if col != 'variant_name']]

    return df

# Process known_variants.tsv to drop specific columns and sort
known_variants_path = os.path.join(output_folder, 'known_variants.tsv')
if os.path.exists(known_variants_path):
    known_variants_df = pd.read_csv(known_variants_path, sep='\t')

    # Drop columns containing '_AD' and others like 'CHROM', 'POS', 'REF', 'ALT'
    columns_to_drop = [col for col in known_variants_df.columns if '_AD' in col or col in ['CHROM', 'POS', 'REF', 'ALT']]
    known_variants_df.drop(columns=columns_to_drop, inplace=True)

    # Rename and sort the remaining columns
    known_variants_df = modify_dataframe(known_variants_df)

    # Save the modified DataFrame
    filtered_output_path = os.path.join(output_folder, 'known_variants_filtered.tsv')
    known_variants_df.to_csv(filtered_output_path, sep='\t', index=False)

    print(f"Filtered file saved to: {filtered_output_path}")
else:
    print(f"Error: {known_variants_path} does not exist.")