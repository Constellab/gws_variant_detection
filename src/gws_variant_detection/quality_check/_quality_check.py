
import sys
import pandas as pd
import os
import subprocess


class NanoPlotQualityChecker:
    def __init__(self, input_folder, output_folder, num_threads):
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.num_threads = num_threads

    def check_quality(self):
        # Crée le dossier de sortie s'il n'existe pas
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        # Liste tous les fichiers FASTQ dans le dossier d'entrée
        fastq_files = [f for f in os.listdir(self.input_folder) if f.endswith('.fastq') or f.endswith('.fastq.gz')]

        # Exécute NanoPlot pour chaque fichier FASTQ
        for fastq_file in fastq_files:
            input_path = os.path.join(self.input_folder, fastq_file)
            output_path = os.path.join(self.output_folder, fastq_file.split('.')[0])  # Construction du chemin de sortie
            nanoplot_cmd = f'NanoPlot --fastq {input_path} -o {output_path} --threads {self.num_threads} --tsv_stats --raw'
            subprocess.run(nanoplot_cmd, shell=True)

class NanoPlotSummary:
    def __init__(self, quality_output_folder):
        self.quality_output_folder = quality_output_folder
        self.desired_columns = [
            'number_of_reads', 'number_of_bases', 'median_read_length', 'mean_read_length',
            'Reads >Q5:', 'Reads >Q7:', 'Reads >Q10:', 'Reads >Q12:', 'Reads >Q15:', 'Reads >Q20:'
        ]

    def clean_column_names(self, df):
        # Supprime les deux-points ":" des noms de colonnes spécifiques
        df.columns = df.columns.str.replace(':', '', regex=True)

    def merge_nano_stats(self):
        # Initialise un DataFrame vide pour stocker les données fusionnées
        merged_data = pd.DataFrame()

        # Recherche récursive des fichiers NanoStats.txt dans le dossier de sortie
        for root, dirs, files in os.walk(self.quality_output_folder):
            for file in files:
                if file == 'NanoStats.txt':
                    # Chemin complet vers le fichier NanoStats.txt
                    nano_stats_file = os.path.join(root, file)

                    # Lecture du fichier NanoStats.txt et fusion des données dans le DataFrame
                    with open(nano_stats_file, 'r') as f:
                        lines = f.readlines()
                        data = {}
                        for line in lines[1:]:
                            parts = line.strip().split('\t')
                            metric = parts[0]
                            value = parts[1].split(')')[0].strip()  # Supprime tout après la parenthèse fermante
                            if metric in self.desired_columns:
                                if metric in [
                                    'Reads >Q5:', 'Reads >Q7:', 'Reads >Q10:', 'Reads >Q12:', 'Reads >Q15:', 'Reads >Q20:'
                                ]:
                                    value += ')'  # Préserve la parenthèse fermante pour ces métriques
                                data[metric] = value
                        # Utilise le nom du dossier comme nom de colonne
                        folder_name = os.path.basename(os.path.dirname(nano_stats_file))
                        # Convertit les données en DataFrame avec une seule colonne
                        if data:
                            merged_data[folder_name] = pd.DataFrame(data.values(), index=data.keys())

        # Transpose le DataFrame pour avoir les codes-barres en lignes et les métriques en colonnes
        merged_data = merged_data.transpose()

        # Affiche les valeurs de l'index pour le débogage
        print("Valeurs de l'index de merged_data :")
        print(merged_data.index)

        # Extrait le numéro du code-barres depuis l'index
        merged_data['Barcode_Number'] = merged_data.index.to_series().str.extract(r'barcode(\d+)')[0]

        # Identifie les lignes avec des numéros de code-barres NaN
        nan_barcodes = merged_data[merged_data['Barcode_Number'].isna()]
        if not nan_barcodes.empty:
            print("Les entrées suivantes ont des numéros de code-barres NaN et seront supprimées :")
            print(nan_barcodes.index.tolist())

        # Supprime les lignes avec des numéros de code-barres NaN
        merged_data = merged_data[merged_data['Barcode_Number'].notna()]

        # Convertit Barcode_Number en entier
        merged_data['Barcode_Number'] = merged_data['Barcode_Number'].astype(int)

        # Trie par numéro de code-barres
        merged_data.sort_values(by='Barcode_Number', inplace=True)
        merged_data.drop(columns=['Barcode_Number'], inplace=True)

        # Nettoie les noms de colonnes
        self.clean_column_names(merged_data)

        # Écrit le DataFrame fusionné dans un seul fichier NanoStats.txt
        output_file = os.path.join('NanoStats_Summary.txt')
        merged_data.to_csv(output_file, sep='\t')



# Example usage
input_folder = sys.argv[1]
num_threads = int(sys.argv[2])
output_folder = sys.argv[3]

# Création de l'instance et exécution de la vérification de qualité
nanoplot_checker = NanoPlotQualityChecker(input_folder, output_folder, num_threads)
nanoplot_checker.check_quality()

# Fusion des fichiers NanoStats.txt en une seule table et tri par numéro de code-barres
summary_generator = NanoPlotSummary(output_folder)
summary_generator.merge_nano_stats()

