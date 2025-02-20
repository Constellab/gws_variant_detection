
import sys
import os
import pandas as pd
import numpy as np


class CooccurrenceMatrixGenerator:
    def __init__(
        self,
        input_csv,
        output_counts_folder,
        output_freq_folder,
        output_subset_folder
    ):

        self.input_csv = input_csv
        self.output_counts_folder = output_counts_folder
        self.output_freq_folder = output_freq_folder
        self.output_subset_folder = output_subset_folder

        self.create_folder(self.output_counts_folder)
        self.create_folder(self.output_freq_folder)
        self.create_folder(self.output_subset_folder)

    def create_folder(self, folder_path):
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

    def load_read_variants(self):

        if not os.path.exists(self.input_csv):
            raise FileNotFoundError(f"Le fichier {self.input_csv} n'existe pas.")

        df = pd.read_csv(self.input_csv)
        if 'Variants' not in df.columns or 'Read' not in df.columns:
            raise ValueError("Le fichier CSV doit contenir 'Read' et 'Variants'.")
        return df

    def extract_variants(self, variants_str):

        try:
            variants_str = variants_str.strip('[]')
            return [
                variant.strip().replace("'", "")
                for variant in variants_str.split(',')
                if variant.strip()
            ]
        except Exception as e:
            print(f"Erreur lors de l'analyse de la chaîne '{variants_str}': {e}")
            return []

    def extract_barcode(self, read_id):

        try:
            return read_id.split('-')[0]
        except Exception as e:
            print(f"Erreur d'extraction du barcode : {e}")
            return "Unknown"


    # (1) GÉNÉRER LES MATRICES 2D (IGNORE 0 OU >=3 VARIANTS)

    def generate_2d_matrices(self, df):

        df['Barcode'] = df['Read'].apply(self.extract_barcode)
        grouped = df.groupby('Barcode')

        for barcode, group in grouped:

            # Rassemble toutes les variantes distinctes
            all_variants = set()
            for variants_str in group['Variants']:
                var_list = self.extract_variants(variants_str)
                all_variants.update(var_list)

            if not all_variants:
                print(f"Aucune variante pour {barcode}. On ignore.")
                continue

            all_variants = sorted(all_variants)
            variant_indices = {v: i for i, v in enumerate(all_variants)}
            count_matrix = np.zeros((len(all_variants), len(all_variants)), dtype=int)
            total_reads = len(group)

            # Parcours des reads
            for variants_str in group['Variants']:
                var_list = self.extract_variants(variants_str)

                if len(var_list) == 1:
                    # Diagonale
                    v = var_list[0]
                    if v in variant_indices:
                        idx = variant_indices[v]
                        count_matrix[idx, idx] += 1
                elif len(var_list) == 2:
                    # Off-diagonal
                    v1, v2 = var_list
                    if v1 in variant_indices and v2 in variant_indices:
                        i1 = variant_indices[v1]
                        i2 = variant_indices[v2]
                        count_matrix[i1, i2] += 1
                        count_matrix[i2, i1] += 1
                else:
                    # 0 ou >=3 => on ignore
                    continue

            # Matrice de comptes
            count_df = pd.DataFrame(count_matrix, index=all_variants, columns=all_variants)
            counts_path = os.path.join(self.output_counts_folder, f"{barcode}_cooccurrence_matrix_counts.csv")
            count_df.to_csv(counts_path)

            # Matrice de fréquences
            freq_matrix = (count_matrix / total_reads) * 100
            freq_matrix = np.round(freq_matrix, 2)
            freq_df = pd.DataFrame(freq_matrix, index=all_variants, columns=all_variants)
            freqs_path = os.path.join(self.output_freq_folder, f"{barcode}_cooccurrence_matrix_frequencies.csv")
            freq_df.to_csv(freqs_path)

    #
    # (2) TABLE DES SUBSETS (y compris 0 variants => 'wildtype')
    #
    def generate_subset_tables(self, df):
        """
        Pour chaque barcode, on produit un CSV listant:
         - Subset_of_variants
         - Count
         - Frequency(%)
        Les reads avec 0 variants => 'wildtype'.
        On retire tous guillemets restants dans les noms de variantes.
        """
        df['Barcode'] = df['Read'].apply(self.extract_barcode)

        counts_dict = {}  # (barcode, subset_tuple) -> count
        total_reads_per_barcode = {}

        for _, row in df.iterrows():
            barcode = row['Barcode']
            var_list = self.extract_variants(row['Variants'])
            # Remove leftover quotes in variant names
            var_list = [v.replace('"', '').replace("'", "") for v in var_list if v.strip()]

            if len(var_list) == 0:
                # => wildtype
                subset_key = ("wildtype",)
            else:
                subset_key = tuple(sorted(var_list))

            if barcode not in total_reads_per_barcode:
                total_reads_per_barcode[barcode] = 0
            total_reads_per_barcode[barcode] += 1

            dict_key = (barcode, subset_key)
            counts_dict[dict_key] = counts_dict.get(dict_key, 0) + 1

        # Un CSV par barcode
        barcodes = sorted(total_reads_per_barcode.keys())
        for bc in barcodes:
            total_bc = total_reads_per_barcode[bc]
            bc_subsets = []
            bc_counts = []
            bc_freqs = []

            for (bcode, subset_tuple), ccount in counts_dict.items():
                if bcode == bc:
                    # Transform subset_tuple en string
                    if len(subset_tuple) == 1 and subset_tuple[0] == "wildtype":
                        subset_str = "wildtype"
                    else:
                        # Join with comma
                        subset_str = ",".join(
                            x.replace('"', '').replace("'", "")
                            for x in subset_tuple
                        )
                    bc_subsets.append(subset_str)
                    bc_counts.append(ccount)
                    freq = round((ccount / total_bc) * 100, 2)
                    bc_freqs.append(freq)

            out_df = pd.DataFrame({
                "Subset_of_variants": bc_subsets,
                "Count": bc_counts,
                "Frequency(%)": bc_freqs
            })
            # On trie par Count desc
            out_df = out_df.sort_values("Count", ascending=False)

            out_path = os.path.join(self.output_subset_folder, f"{bc}_subset_cooccurrence.csv")
            out_df.to_csv(out_path, index=False)

    def run(self):
        df = self.load_read_variants()
        self.generate_2d_matrices(df)
        self.generate_subset_tables(df)

input_csv = sys.argv[1]
output_counts_folder = sys.argv[2]
output_freq_folder = sys.argv[3]
output_subset_folder = sys.argv[4]

generator = CooccurrenceMatrixGenerator(
    input_csv=input_csv,
    output_counts_folder=output_counts_folder,
    output_freq_folder=output_freq_folder,
    output_subset_folder=output_subset_folder
)
generator.run()