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
        subprocess.run(['bcftools', 'norm', '-f', self.ref_genome, '-o', final_vcf, intermediate_vcf], check=True)

        self.extract_data_to_tsv(final_vcf, vcf_file)

    def extract_data_to_tsv(self, final_vcf_path, vcf_file):
        output_tsv = os.path.join(self.output_folder, f"{os.path.splitext(vcf_file)[0]}_variants.tsv")
        with open(final_vcf_path, 'r') as vcf, open(output_tsv, 'w') as tsv:
            tsv.write("CHROM\tPOS\tREF\tALT\tAD\tVARIANT_FREQ\tVARIANT_COUNT\tQUAL\n")
            for line in vcf:
                if not line.startswith('#'):
                    columns = line.strip().split('\t')
                    qual = float(columns[5])
                    if qual >= self.min_qual:
                        format_parts = columns[8].split(':')
                        ad_index = format_parts.index('AD') if 'AD' in format_parts else None
                        if ad_index is not None and len(columns) > 9:
                            sample_info = columns[9].split(':')
                            ad_value = sample_info[ad_index]
                            ad_parts = ad_value.split(',')
                            if len(ad_parts) == 2 and all(x.isdigit() for x in ad_parts):
                                ref_count, alt_count = map(int, ad_parts)
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
        merged_df = merged_df.loc[:,~merged_df.columns.duplicated()]
        merged_df = merged_df.groupby(['CHROM', 'POS', 'REF', 'ALT']).first().reset_index()
        merged_df.to_csv(os.path.join(self.output_folder, 'merged_variants.tsv'), sep='\t', index=False)
        return merged_df

    def parse_variants(self, merged_df):
        # Define known variants for each gene
        gene_variants = {
            'Cytb': [
                {'variant_name': 'Variant4_S34L', 'CHROM': 'Cytb', 'POS': 919, 'REF': 'C', 'ALT': 'T'},
                {'variant_name': 'Variant5_F129L_A', 'CHROM': 'Cytb', 'POS': 1205, 'REF': 'T', 'ALT': 'A'},
                {'variant_name': 'Variant5_F129L_G', 'CHROM': 'Cytb', 'POS': 1205, 'REF': 'T', 'ALT': 'G'},
                {'variant_name': 'Variant6_G143A', 'CHROM': 'Cytb', 'POS': 1246, 'REF': 'G', 'ALT': 'C'},
                {'variant_name': 'Variant3_L201S', 'CHROM': 'Cytb', 'POS': 1420, 'REF': 'T', 'ALT': 'C'},
                {'variant_name': 'Variant1_E203-DE-V204', 'CHROM': 'Cytb', 'POS': 1422, 'REF': 'C', 'ALT': 'CATGAGG'},
                {'variant_name': 'Variant2_E203-VE-V204', 'CHROM': 'Cytb', 'POS': 1431, 'REF': 'G', 'ALT': 'GAGGTAG'}
            ],
            'OSBP': [
                {'variant_name': 'G770V', 'CHROM': 'OSBP', 'POS': 3889, 'REF': 'G', 'ALT': 'T'},
                {'variant_name': 'N837I', 'CHROM': 'OSBP', 'POS': 4090, 'REF': 'A', 'ALT': 'T'},
                {'variant_name': 'N837I_ATT', 'CHROM': 'OSBP', 'POS': 4090, 'REF': 'A', 'ALT': 'ATT'},
                {'variant_name': 'N837I_ATC', 'CHROM': 'OSBP', 'POS': 4090, 'REF': 'A', 'ALT': 'ATC'},
                {'variant_name': 'L863W', 'CHROM': 'OSBP', 'POS': 4168, 'REF': 'T', 'ALT': 'G'}
            ],
            'Tubulin': [
                {'variant_name': 'C239S_AGC', 'CHROM': 'Tubulin', 'POS': 1541, 'REF': 'T', 'ALT': 'A'},
                {'variant_name': 'C239G_GGC', 'CHROM': 'Tubulin', 'POS': 1541, 'REF': 'T', 'ALT': 'G'},
                {'variant_name': 'C239S_TCC', 'CHROM': 'Tubulin', 'POS': 1542, 'REF': 'G', 'ALT': 'C'}
            ],
            'ATPase': []
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

        # Reorder columns to place variant_name first
        cols = ['variant_name'] + [col for col in merged_known.columns if col != 'variant_name']
        merged_known = merged_known[cols]


        merged_known.to_csv(os.path.join(self.output_folder, 'known_variants.tsv'), sep='\t', index=False)
        unknown = merged_df[~merged_df.set_index(['CHROM', 'POS', 'REF', 'ALT']).index.isin(known_df.set_index(['CHROM', 'POS', 'REF', 'ALT']).index)]
        for col in unknown.columns:
            if '_QUAL' in col:
                unknown.drop(col, axis=1, inplace=True)

        unknown.to_csv(os.path.join(self.output_folder, 'unknown_variants.tsv'), sep='\t', index=False)

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


