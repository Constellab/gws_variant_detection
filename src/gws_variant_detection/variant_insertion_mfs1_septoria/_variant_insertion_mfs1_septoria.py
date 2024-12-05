import sys
import os
import pysam
import mappy as mp
import pandas as pd
import re

class BAMClassifier:
    def __init__(self, bam_dir, output_dir, combined_reference_fasta, identity_threshold, coverage_threshold):
        self.bam_dir = bam_dir
        self.output_dir = output_dir
        self.combined_reference_fasta = combined_reference_fasta
        self.identity_threshold = identity_threshold
        self.coverage_threshold = coverage_threshold

    def classify_samples(self):
        # Lister les fichiers BAM dans le répertoire
        bam_files = [f for f in os.listdir(self.bam_dir) if f.endswith(".bam")]
        summary_data = []

        # S'assurer que les répertoires de sortie existent
        os.makedirs(self.output_dir, exist_ok=True)

        # Charger les séquences de référence combinées et créer un aligner
        aligner = mp.Aligner(self.combined_reference_fasta, preset="map-ont")
        if not aligner:
            raise ValueError("Échec du chargement des séquences de référence combinées.")

        for bam_file in bam_files:
            bam_path = os.path.join(self.bam_dir, bam_file)
            sample_name = os.path.splitext(bam_file)[0]

            print(f"Traitement de {bam_path}")

            # Classifier les lectures dans le fichier BAM
            classified_reads, total_reads = self._classify_reads(bam_path, aligner)

            # Écrire les lectures classées dans de nouveaux fichiers BAM
            for ref_name, reads in classified_reads.items():
                if reads:  # Écrire seulement s'il y a des lectures
                    bam_output_path = os.path.join(self.output_dir, ref_name, f"{sample_name}_{ref_name}.bam")
                    self._write_bam(reads, bam_path, bam_output_path)

            # Préparer les données pour le résumé
            data = {'Sample': sample_name, 'Total_Reads': total_reads}
            for ref_name in aligner.seq_names:
                count = len(classified_reads.get(ref_name, []))
                frequency = (count / total_reads * 100) if total_reads > 0 else 0
                data[f"{ref_name}_Count"] = count
                data[f"{ref_name}_Frequency"] = round(frequency, 2)
            # Compte pour la catégorie 'Other'
            other_count = len(classified_reads.get('Other', []))
            other_frequency = (other_count / total_reads * 100) if total_reads > 0 else 0
            data['Other_Count'] = other_count
            data['Other_Frequency'] = round(other_frequency, 2)

            summary_data.append(data)

        # Sauvegarder les données du résumé dans un fichier CSV
        summary_df = pd.DataFrame(summary_data)

        # Fonction pour extraire le numéro de barcode du nom de l'échantillon
        def extract_barcode_number(sample_name):
            match = re.search(r'barcode(\d+)', sample_name)
            if match:
                return int(match.group(1))  # Extraire la partie numérique du barcode en tant qu'entier
            return float('inf')  # Si aucun barcode n'est trouvé, le placer à la fin lors du tri

        # Trier le DataFrame par numéro de barcode extrait
        summary_df['Barcode_Num'] = summary_df['Sample'].apply(extract_barcode_number)
        summary_df = summary_df.sort_values('Barcode_Num').drop('Barcode_Num', axis=1)

        # Écrire le DataFrame trié dans un fichier CSV
        summary_csv_path = os.path.join(self.output_dir, "classification_summary.csv")
        summary_df.to_csv(summary_csv_path, index=False)
        print(f"Résumé sauvegardé dans {summary_csv_path}")

    def _classify_reads(self, bam_path, aligner):
        classified_reads = {ref_name: [] for ref_name in aligner.seq_names}
        classified_reads['Other'] = []
        total_reads = 0

        with pysam.AlignmentFile(bam_path, "rb") as bam_file:
            for read in bam_file:
                if read.is_secondary or read.is_supplementary:
                    continue

                total_reads += 1

                seq = read.query_sequence

                if seq is None:
                    continue

                # Aligner la séquence de la lecture sur les références combinées
                hits = list(aligner.map(seq))

                if not hits:
                    classified_reads['Other'].append(read)
                    continue

                # Trier les hits par score de mapping (mapq)
                hits.sort(key=lambda x: x.mapq, reverse=True)
                best_hit = hits[0]

                # Calculer le pourcentage d'identité
                identity = (best_hit.mlen / best_hit.blen) * 100

                # Obtenir la longueur de la séquence de référence correspondante
                ref_seq = aligner.seq(best_hit.ctg)
                ref_length = len(ref_seq)

                # Calculer la couverture alignée de la référence
                ref_coverage = ((best_hit.r_en - best_hit.r_st) / ref_length) * 100

                # Vérifier les critères de seuil
                if identity >= self.identity_threshold and ref_coverage >= self.coverage_threshold:
                    ref_name = best_hit.ctg
                    classified_reads[ref_name].append(read)
                else:
                    classified_reads['Other'].append(read)

        return classified_reads, total_reads

    def _write_bam(self, reads, original_bam_path, bam_output_path):
        # Écrire les lectures classées dans un nouveau fichier BAM
        os.makedirs(os.path.dirname(bam_output_path), exist_ok=True)  # S'assurer que le répertoire existe
        with pysam.AlignmentFile(original_bam_path, "rb") as original_bam:
            with pysam.AlignmentFile(bam_output_path, "wb", header=original_bam.header) as out_bam:
                for read in reads:
                    out_bam.write(read)

        # Trier et indexer le fichier BAM
        sorted_bam_path = bam_output_path.replace(".bam", ".sorted.bam")
        pysam.sort("-o", sorted_bam_path, bam_output_path)
        pysam.index(sorted_bam_path)

        # Supprimer le fichier BAM non trié
        os.remove(bam_output_path)

    # Définir les répertoires d'entrée et de sortie, et le chemin vers la séquence de référence combinée
bam_dir = sys.argv[1]
output_dir = sys.argv[2]
combined_reference_fasta = sys.argv[3]
identity_threshold = int(sys.argv[4])
coverage_threshold = int(sys.argv[5])

classifier = BAMClassifier(
    bam_dir,
    output_dir,
    combined_reference_fasta,
    identity_threshold,
    coverage_threshold
)
classifier.classify_samples()
