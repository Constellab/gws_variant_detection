import sys
import os
import subprocess
import gzip
import mappy as mp
import pandas as pd

class FastpFiltering:
    def __init__(self, input_folder, fastp_output_folder, identity_coverage_output_folder, final_output_folder, summary_output_folder,
                 num_threads, qualified_quality_phred, reads_to_process, min_identity, min_coverage,
                 reference_genome, skip_read_threshold, min_read_length):
        self.input_folder = input_folder
        self.fastp_output_folder = fastp_output_folder
        self.identity_coverage_output_folder = identity_coverage_output_folder
        self.final_output_folder = final_output_folder
        self.summary_output_folder = summary_output_folder
        self.num_threads = num_threads
        self.qualified_quality_phred = qualified_quality_phred
        self.reads_to_process = reads_to_process  # Nombre final de reads à sélectionner (5000)
        self.min_identity = min_identity
        self.min_coverage = min_coverage
        self.reference_genome = reference_genome
        self.skip_read_threshold = skip_read_threshold
        self.min_read_length = min_read_length  # Longueur minimale de read (400 pb)
        self.aligner = mp.Aligner(reference_genome, preset="map-ont")

        # Création des dossiers de sortie si nécessaire
        os.makedirs(self.fastp_output_folder, exist_ok=True)
        os.makedirs(self.identity_coverage_output_folder, exist_ok=True)
        os.makedirs(self.final_output_folder, exist_ok=True)
        os.makedirs(self.summary_output_folder, exist_ok=True)

    def count_reads_in_fastq(self, fastq_file):
        """Compte le nombre de reads dans un fichier FASTQ."""
        try:
            cmd = f'zcat {fastq_file} | wc -l'
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
            num_lines = int(result.stdout.strip())
            return num_lines // 4  # Chaque read occupe 4 lignes dans un fichier FASTQ
        except subprocess.CalledProcessError as e:
            print(f"Error counting reads in {fastq_file}: {e}", file=sys.stderr)
            return 0

    def check_alignment_coverage(self, read_sequence):
        """Vérifie qu'un read respecte le seuil d'identité et de couverture lors de l'alignement sur la référence."""
        for hit in self.aligner.map(read_sequence):
            query_length = hit.q_en - hit.q_st  # Longueur du read

            if query_length == 0 or hit.blen == 0:  # Éviter une division par zéro
                continue

            identity = (hit.mlen / hit.blen) * 100  # Pourcentage de bases identiques dans l'alignement
            coverage = (hit.blen / query_length) * 100  # Pourcentage du read aligné

            if identity >= self.min_identity and coverage >= self.min_coverage:
                return identity, coverage, query_length  # Read validé

        return None, None, None  # Read rejeté

    def filter_reads(self):
        """Filtre les reads selon la qualité fastp, la longueur minimale et les seuils d'identité/couverture.
           Enfin, limite le nombre de reads finaux à 'reads_to_process' (ici 5000)."""
        summary_data = []

        fastq_files = [f for f in os.listdir(self.input_folder) if f.endswith('.fastq.gz')]

        for fastq_file in fastq_files:
            input_path = os.path.join(self.input_folder, fastq_file)
            fastp_output_path = os.path.join(self.fastp_output_folder, fastq_file)
            identity_coverage_output_path = os.path.join(self.identity_coverage_output_folder, fastq_file)
            final_output_path = os.path.join(self.final_output_folder, fastq_file)

            # Étape 1 : Exécution de fastp pour le filtrage qualité
            fastp_cmd = f'gunzip -c {input_path} | chopper -q {self.qualified_quality_phred} --threads {self.num_threads} | gzip > {fastp_output_path}'
            result = subprocess.run(fastp_cmd, shell=True)
            if result.returncode != 0:
                print(f"Error with fastp on {fastq_file}", file=sys.stderr)
                continue

            # Étape 2 : Filtrage des reads basé sur la qualité, la longueur minimale, et l'identité/couverture
            valid_reads = []
            final_reads = []

            with gzip.open(fastp_output_path, 'rt') as infile:
                fastq_lines = infile.readlines()

            with gzip.open(identity_coverage_output_path, 'wt') as outfile:
                for i in range(0, len(fastq_lines), 4):
                    read_name = fastq_lines[i].strip()
                    read_sequence = fastq_lines[i + 1].strip()
                    read_quality = fastq_lines[i + 3].strip()

                    # Vérifier la longueur minimale du read
                    if len(read_sequence) < self.min_read_length:
                        continue

                    # Le read passe le filtre fastp et la condition de longueur
                    valid_reads.append((read_name, read_sequence, read_quality))

                    identity, coverage, read_length = self.check_alignment_coverage(read_sequence)
                    if identity is not None and coverage is not None:
                        outfile.write(f"{read_name}\n{read_sequence}\n+\n{read_quality}\n")
                        final_reads.append((read_name, read_sequence, read_quality, identity, coverage, read_length))

            # À ce stade, 'final_reads' contient tous les reads qui passent les filtres de qualité, longueur et alignement.
            # Appliquer ensuite la limite de 5000 reads, c'est-à-dire le critère doit être appliqué en dernier.
            if len(final_reads) > self.reads_to_process:
                # On prend les 5000 premiers reads
                final_reads = final_reads[:self.reads_to_process]

            # Écriture du fichier final de sortie
            with gzip.open(final_output_path, 'wt') as final_file:
                for read in final_reads:
                    final_file.write(f"{read[0]}\n{read[1]}\n+\n{read[2]}\n")

            # Informations récapitulatives
            summary_data.append({
                "File": fastq_file,
                "Total_Reads_After_Phred_and_Length": len(valid_reads),
                "Total_Reads_After_Identity_Coverage": len(final_reads),
                "Final_Filtered_Reads": len(final_reads),
                "Min_Identity": self.min_identity,
                "Min_Coverage": self.min_coverage,
                "Phred_Threshold": self.qualified_quality_phred,
                "Min_Read_Length": self.min_read_length,
                "Reads_Selected": self.reads_to_process if len(final_reads) >= self.reads_to_process else len(final_reads)
            })

        # Sauvegarde du résumé en CSV
        summary_df = pd.DataFrame(summary_data)
        summary_csv_path = os.path.join(self.summary_output_folder, "filtering_summary.csv")
        summary_df.to_csv(summary_csv_path, index=False)
        print(f"Summary saved to {summary_csv_path}")


# Paramètres de la ligne de commande (exemple d'utilisation)
input_folder = sys.argv[1]
fastp_output_folder = sys.argv[2]
identity_coverage_output_folder = sys.argv[3]
final_output_folder = sys.argv[4]
summary_output_folder = sys.argv[5]

num_threads = int(sys.argv[6])
qualified_quality_phred = int(sys.argv[7])
reads_to_process = int(sys.argv[8])
min_identity = int(sys.argv[9])
min_coverage = int(sys.argv[10])
reference_genome = sys.argv[11]
skip_read_threshold = int(sys.argv[12])
min_read_length = int(sys.argv[13])  

fastp_filtering = FastpFiltering(
    input_folder, fastp_output_folder, identity_coverage_output_folder, final_output_folder, summary_output_folder,
    num_threads, qualified_quality_phred, reads_to_process, min_identity, min_coverage, reference_genome, skip_read_threshold, min_read_length
)
fastp_filtering.filter_reads()
