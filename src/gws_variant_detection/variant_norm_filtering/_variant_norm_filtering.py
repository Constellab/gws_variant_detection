import os
import pandas as pd
import subprocess
import sys

class VCFNormalizer:
    def __init__(self, vcf_folder, ref_genome, gene_name, reads_to_process, final_normalized_output_folder, output_folder, min_qual):
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
        subprocess.run(['bcftools', 'norm', '-f', self.ref_genome, '-cw', '-o', final_vcf, intermediate_vcf], check=True)

        self.extract_data_to_tsv(final_vcf, vcf_file)

    def extract_data_to_tsv(self, final_vcf_path, vcf_file):
        output_tsv = os.path.join(self.output_folder, f"{os.path.splitext(vcf_file)[0]}_variants.tsv")

        if self.gene_name == 'Tubulin':
            self.extract_data_tubulin(final_vcf_path, vcf_file, output_tsv)
        else:
            self.extract_data_cytb_osbp(final_vcf_path, vcf_file, output_tsv)

    def extract_data_tubulin(self, final_vcf_path, vcf_file, output_tsv):
        # Variables to store counts and frequencies for specific positions
        variant_data = {
            '1541_A': {'ref_count': 0, 'alt_count': 0, 'variant_freq': 0},
            '1541_G': {'ref_count': 0, 'alt_count': 0, 'variant_freq': 0}
        }

        # First pass to extract REF and ALT counts for specific positions
        with open(final_vcf_path, 'r') as vcf:
            for line in vcf:
                if not line.startswith('#'):
                    columns = line.strip().split('\t')
                    pos = int(columns[1])
                    ref = columns[3]
                    alt = columns[4]
                    format_parts = columns[8].split(':')
                    ad_index = format_parts.index('AD') if 'AD' in format_parts else None
                    if ad_index is not None and len(columns) > 9:
                        sample_info = columns[9].split(':')
                        ad_value = sample_info[ad_index]
                        ad_parts = ad_value.split(',')
                        if len(ad_parts) == 2 and all(x.isdigit() for x in ad_parts):
                            if pos == 1541 and ref == 'T' and alt == 'A':
                                variant_data['1541_A']['ref_count'] = int(ad_parts[0])
                                variant_data['1541_A']['alt_count'] = int(ad_parts[1])
                                total_count_1541_A = variant_data['1541_A']['ref_count'] + variant_data['1541_A']['alt_count']
                                variant_data['1541_A']['variant_freq'] = (variant_data['1541_A']['alt_count'] / total_count_1541_A * 100) if total_count_1541_A > 0 else 0
                            elif pos == 1541 and ref == 'T' and alt == 'G':
                                variant_data['1541_G']['ref_count'] = int(ad_parts[0])
                                variant_data['1541_G']['alt_count'] = int(ad_parts[1])
                                total_count_1541_G = variant_data['1541_G']['ref_count'] + variant_data['1541_G']['alt_count']
                                variant_data['1541_G']['variant_freq'] = (variant_data['1541_G']['alt_count'] / total_count_1541_G * 100) if total_count_1541_G > 0 else 0

        # If the gene is Tubulin, normalize frequencies if combined exceeds 100%
        if variant_data['1541_A']['variant_freq'] + variant_data['1541_G']['variant_freq'] > 100:
            total_ref_alt_1541 = variant_data['1541_A']['ref_count'] + variant_data['1541_A']['alt_count'] + variant_data['1541_G']['ref_count'] + variant_data['1541_G']['alt_count']
            if total_ref_alt_1541 > 0:
                ratio_1541_A = (variant_data['1541_A']['ref_count'] + variant_data['1541_A']['alt_count']) / total_ref_alt_1541
                ratio_1541_G = (variant_data['1541_G']['ref_count'] + variant_data['1541_G']['alt_count']) / total_ref_alt_1541
                print(f"File: {vcf_file}")
                print(f"Ratio for 1541_A: {ratio_1541_A}, FREQ for 1541_A: {variant_data['1541_A']['variant_freq']}")
                print(f"Ratio for 1541_G: {ratio_1541_G}, FREQ for 1541_G: {variant_data['1541_G']['variant_freq']}")
                variant_data['1541_A']['variant_freq'] = round(ratio_1541_A * 100, 2)
                variant_data['1541_G']['variant_freq'] = round(ratio_1541_G * 100, 2)
                print(f"Normalized FREQ for 1541_A: {variant_data['1541_A']['variant_freq']}")
                print(f"Normalized FREQ for 1541_G: {variant_data['1541_G']['variant_freq']}")

        # Second pass to extract other data
        with open(final_vcf_path, 'r') as vcf, open(output_tsv, 'w') as tsv:
            tsv.write("CHROM\tPOS\tREF\tALT\tAD\tVARIANT_FREQ\tVARIANT_COUNT\tQUAL\n")
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
                                variant_freq = 0  # Initialize variant_freq to avoid UnboundLocalError
                                if pos == 1541 and ref_count > 0 and alt_count > 0:
                                    if columns[3] == 'T' and columns[4] == 'A':
                                        variant_freq = variant_data['1541_A']['variant_freq']
                                    elif columns[3] == 'T' and columns[4] == 'G':
                                        variant_freq = variant_data['1541_G']['variant_freq']
                                    variant_count = round((variant_freq / 100) * self.reads_to_process)
                                    tsv.write(f"{columns[0]}\t{pos}\t{columns[3]}\t{columns[4]}\t{ad_value}\t{variant_freq:.2f}%\t{variant_count}\t{qual}\n")
                                elif pos == 1542:
                                    variant_freq = (alt_count / total_count * 100) if total_count > 0 else 0
                                    variant_count = round((variant_freq / 100) * self.reads_to_process)
                                    tsv.write(f"{columns[0]}\t{pos}\t{columns[3]}\t{columns[4]}\t{ad_value}\t{variant_freq:.2f}%\t{variant_count}\t{qual}\n")

    def extract_data_cytb_osbp(self, final_vcf_path, vcf_file, output_tsv):
        ref_count_1422 = alt_count_1422 = 0
        ref_count_1428 = alt_count_1428 = 0
        ref_count_1427 = alt_count_1427 = 0
        ref_count_1430 = alt_count_1430 = 0

        # First pass to extract REF and ALT counts for specific positions
        with open(final_vcf_path, 'r') as vcf:
            for line in vcf:
                if not line.startswith('#'):
                    columns = line.strip().split('\t')
                    pos = int(columns[1])
                    ref = columns[3]
                    alt = columns[4]
                    format_parts = columns[8].split(':')
                    ad_index = format_parts.index('AD') if 'AD' in format_parts else None
                    if ad_index is not None and len(columns) > 9:
                        sample_info = columns[9].split(':')
                        ad_value = sample_info[ad_index]
                        ad_parts = ad_value.split(',')
                        if len(ad_parts) == 2 and all(x.isdigit() for x in ad_parts):
                            if pos == 1422 and ref == 'C' and alt == 'CATGAGG':
                                ref_count_1422 = int(ad_parts[0])
                                alt_count_1422 = int(ad_parts[1])
                            elif pos == 1428 and ref == 'G' and alt == 'GA':
                                ref_count_1428 = int(ad_parts[0])
                                alt_count_1428 = int(ad_parts[1])
                            elif pos == 1427 and ref == 'G' and alt == 'A':
                                ref_count_1427 = int(ad_parts[0])
                                alt_count_1427 = int(ad_parts[1])
                            elif pos == 1430 and ref == 'A' and alt == 'G':
                                ref_count_1430 = int(ad_parts[0])
                                alt_count_1430 = int(ad_parts[1])

        # Calculate total counts for each position
        total_count_1422 = ref_count_1422 + alt_count_1422
        total_count_1428 = ref_count_1428 + alt_count_1428
        total_count_1427 = ref_count_1427 + alt_count_1427
        total_count_1430 = ref_count_1430 + alt_count_1430

        # Calculate variant frequencies using the new method
        common_denominator = total_count_1422 * total_count_1428 * total_count_1427 * total_count_1430

        if common_denominator > 0:
            variant_freq_1422 = (alt_count_1422 * total_count_1428 * total_count_1427 * total_count_1430) / common_denominator * 100
            variant_freq_1428 = (alt_count_1428 * total_count_1422 * total_count_1427 * total_count_1430) / common_denominator * 100
            variant_freq_1427 = (alt_count_1427 * total_count_1422 * total_count_1428 * total_count_1430) / common_denominator * 100
            variant_freq_1430 = (alt_count_1430 * total_count_1422 * total_count_1428 * total_count_1427) / common_denominator * 100
        else:
            variant_freq_1422 = variant_freq_1428 = variant_freq_1427 = variant_freq_1430 = 0


        # Second pass to extract other data and process POS 1422 specially
        with open(final_vcf_path, 'r') as vcf, open(output_tsv, 'w') as tsv:
            tsv.write("CHROM\tPOS\tREF\tALT\tAD\tVARIANT_FREQ\tVARIANT_COUNT\tQUAL\n")
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

                                if pos == 1422:
                                    variant_freq = variant_freq_1422 + variant_freq_1428 + variant_freq_1427 + variant_freq_1430
                                    variant_count = round((variant_freq / 100) * self.reads_to_process)
                                else:
                                    total_count = ref_count + alt_count
                                    variant_freq = (alt_count / total_count * 100) if total_count > 0 else 0
                                    variant_count = round((variant_freq / 100) * self.reads_to_process)

                                selected_data = [columns[i] for i in [0, 1, 3, 4]] + [ad_value, f"{variant_freq:.2f}%", f"{variant_count}", str(qual)]
                                tsv.write("\t".join(selected_data) + '\n')

    def merge_variant_files(self):
        tsv_files = [f for f in os.listdir(self.output_folder) if f.endswith('_variants.tsv')]
        dfs = []
        for tsv in tsv_files:
            sample_name = tsv.split('-')[0]
            df = pd.read_csv(os.path.join(self.output_folder, tsv), sep='\t')
            df.columns = [f"{sample_name}_{col}" if col not in ['CHROM', 'POS', 'REF', 'ALT'] else col for col in df.columns]
            dfs.append(df)
        merged_df = dfs[0]
        for df in dfs[1:]:
            merged_df = pd.merge(merged_df, df, on=['CHROM', 'POS', 'REF', 'ALT'], how='outer')
        merged_df.fillna(0, inplace=True)
        merged_df = merged_df.loc[:, ~merged_df.columns.duplicated()]
        merged_df = merged_df.groupby(['CHROM', 'POS', 'REF', 'ALT']).first().reset_index()
        merged_df.to_csv(os.path.join(self.output_folder, 'merged_variants.tsv'), sep='\t', index=False)
        return merged_df

    def parse_variants(self, merged_df):
        # Define known variants for each gene
        gene_variants = {
            'Cytb': [
                {'variant_name': 'S34L', 'CHROM': 'Cytb', 'POS': 919, 'REF': 'C', 'ALT': 'T'},
                {'variant_name': 'F129L_A', 'CHROM': 'Cytb', 'POS': 1205, 'REF': 'T', 'ALT': 'A'},
                {'variant_name': 'F129L_G', 'CHROM': 'Cytb', 'POS': 1205, 'REF': 'T', 'ALT': 'G'},
                {'variant_name': 'G143A', 'CHROM': 'Cytb', 'POS': 1246, 'REF': 'G', 'ALT': 'C'},
                {'variant_name': 'L201S', 'CHROM': 'Cytb', 'POS': 1420, 'REF': 'T', 'ALT': 'C'},
                {'variant_name': 'E203-DE-V204', 'CHROM': 'Cytb', 'POS': 1422, 'REF': 'C', 'ALT': 'CATGAGG'},
                {'variant_name': 'E203-VE-V204', 'CHROM': 'Cytb', 'POS': 1424, 'REF': 'T', 'ALT': 'TGAGGTA'}
            ],
            'OSBP': [
                {'variant_name': 'L863W', 'CHROM': 'OSBP', 'POS': 4169, 'REF': 'T', 'ALT': 'G'},
                {'variant_name': 'N837I', 'CHROM': 'OSBP', 'POS': 4091, 'REF': 'A', 'ALT': 'T'},
                {'variant_name': 'G770V', 'CHROM': 'OSBP', 'POS': 3890, 'REF': 'G', 'ALT': 'T'}
            ],
            'Tubulin': [
                {'variant_name': 'C239S_AGC', 'CHROM': 'Tubulin', 'POS': 1541, 'REF': 'T', 'ALT': 'A'},
                {'variant_name': 'C239G_GGC', 'CHROM': 'Tubulin', 'POS': 1541, 'REF': 'T', 'ALT': 'G'},
                {'variant_name': 'C239S_TCC', 'CHROM': 'Tubulin', 'POS': 1542, 'REF': 'G', 'ALT': 'C'}
            ],
            'ATPase': [
                {'variant_name': 'G759S', 'CHROM': 'ATPase', 'POS': 2857, 'REF': 'G', 'ALT': 'A'},
                {'variant_name': 'I751V', 'CHROM': 'ATPase', 'POS': 2833, 'REF': 'A', 'ALT': 'G'},
                {'variant_name': 'G750A', 'CHROM': 'ATPase', 'POS': 2831, 'REF': 'G', 'ALT': 'C'},
                {'variant_name': 'G750D', 'CHROM': 'ATPase', 'POS': 2831, 'REF': 'G', 'ALT': 'A'},
                {'variant_name': 'G750S', 'CHROM': 'ATPase', 'POS': 2830, 'REF': 'G', 'ALT': 'A'}

            ]
        }

        # Select the known variants for the gene specified at runtime
        known_variants = gene_variants.get(self.gene_name, [])
        known_df = pd.DataFrame(known_variants)
        merged_known = pd.merge(known_df, merged_df, on=['CHROM', 'POS', 'REF', 'ALT'], how='left')

        for col in merged_known.columns:
            if '_AD' in col or '_VARIANT_FREQ' in col or '_VARIANT_COUNT' in col or '_QUAL' in col:
                if 'FREQ' in col:
                    merged_known[col] = merged_known[col].fillna('0%')
                else:
                    merged_known[col] = merged_known[col].fillna(0)

        for col in merged_known.columns:
            if '_QUAL' in col:
                merged_known.drop(col, axis=1, inplace=True)

        # Normalize the frequencies for Tubulin if combined exceeds 100%
        if self.gene_name == 'Tubulin':
            tubulin_variants = ['C239S_AGC', 'C239G_GGC']
            for variant in tubulin_variants:
                freq_cols = [col for col in merged_known.columns if '_VARIANT_FREQ' in col]
                for col in freq_cols:
                    combined_freq = merged_known.loc[merged_known['variant_name'].isin(tubulin_variants), col].apply(lambda x: float(str(x).strip('%'))).sum()
                    if combined_freq > 100:
                        barcode = col.split('_')[0]
                        total_ref_alt_1541 = (
                            merged_known.loc[merged_known['variant_name'] == 'C239S_AGC', f'{barcode}_VARIANT_COUNT'].sum() +
                            merged_known.loc[merged_known['variant_name'] == 'C239G_GGC', f'{barcode}_VARIANT_COUNT'].sum()
                        )
                        if total_ref_alt_1541 > 0:
                            ratio_1541_A = merged_known.loc[merged_known['variant_name'] == 'C239S_AGC', f'{barcode}_VARIANT_COUNT'].sum() / total_ref_alt_1541
                            ratio_1541_G = merged_known.loc[merged_known['variant_name'] == 'C239G_GGC', f'{barcode}_VARIANT_COUNT'].sum() / total_ref_alt_1541
                            if 'C239S_AGC' in merged_known['variant_name'].values:
                                normalized_freq = round(ratio_1541_A * 100, 2)
                                merged_known.loc[merged_known['variant_name'] == 'C239S_AGC', col] = f"{normalized_freq:.2f}%"
                            if 'C239G_GGC' in merged_known['variant_name'].values:
                                normalized_freq = round(ratio_1541_G * 100, 2)
                                merged_known.loc[merged_known['variant_name'] == 'C239G_GGC', col] = f"{normalized_freq:.2f}%"

        # Reorder columns to place variant_name first
        cols = ['variant_name'] + [col for col in merged_known.columns if col != 'variant_name']
        merged_known = merged_known[cols]

        merged_known.to_csv(os.path.join(self.output_folder, 'known_variants.tsv'), sep='\t', index=False)

        unknown = merged_df[~merged_df.set_index(['CHROM', 'POS', 'REF', 'ALT']).index.isin(known_df.set_index(['CHROM', 'POS', 'REF', 'ALT']).index)]
        for col in unknown.columns:
            if '_QUAL' in col:
                unknown.drop(col, axis=1, inplace=True)

        unknown.to_csv(os.path.join(self.output_folder, 'unknown_variants.tsv'), sep='\t', index=False)


vcf_folder = sys.argv[1]
ref_genome = sys.argv[2]
gene_name = str(sys.argv[3])  # Modify as needed for the gene of interest
reads_to_process = int(sys.argv[4])
final_normalized_output_folder = sys.argv[5]
output_folder = sys.argv[6]
min_qual = 0


normalizer = VCFNormalizer(vcf_folder, ref_genome, gene_name, reads_to_process, final_normalized_output_folder, output_folder, min_qual)
normalizer.normalize_vcf_files()
merged_df = normalizer.merge_variant_files()
normalizer.parse_variants(merged_df)

# Cleanup
for file in os.listdir(output_folder):
    if file.endswith(".vcf"):
        os.remove(os.path.join(output_folder, file))

# Helper functions to sort columns
def extract_method_from_column(column_name):
    parts = column_name.split('_')
    try:
        barcode_number = int(parts[0].replace('barcode', ''))
    except ValueError:
        return (float('inf'), '')  # Ensure non-standard columns are sorted last
    method = parts[-1]
    return barcode_number, method

def modify_dataframe(df):
    # Exclude 'variant_name' column before sorting
    other_columns = [col for col in df.columns if col != 'variant_name']
    # Sort columns by barcode number and method
    sorted_columns = sorted(other_columns, key=lambda x: extract_method_from_column(x))
    # Reorder DataFrame columns
    df = df[['variant_name'] + sorted_columns]
    return df

# Process known_variants.tsv to delete specific columns and sort
known_variants_path = os.path.join(output_folder, 'known_variants.tsv')
if os.path.exists(known_variants_path):
    known_variants_df = pd.read_csv(known_variants_path, sep='\t')
    # Drop columns that contain _AD in their names
    columns_to_drop = [col for col in known_variants_df.columns if '_AD' in col]
    known_variants_df.drop(columns=columns_to_drop, inplace=True)
    # Drop columns CHROM, POS, REF, ALT
    known_variants_df.drop(columns=['CHROM', 'POS', 'REF', 'ALT'], inplace=True)
    # Sort columns by barcode number and method
    known_variants_df = modify_dataframe(known_variants_df)
    # Save the modified dataframe
    known_variants_df.to_csv(os.path.join(output_folder, 'known_variants_filtered.tsv'), sep='\t', index=False)


