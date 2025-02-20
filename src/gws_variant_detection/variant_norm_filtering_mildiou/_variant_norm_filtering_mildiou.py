import os
import pandas as pd
import subprocess
import numpy as np
import re
import pysam
import shutil
from multiprocessing import Pool
import ast
import sys

class VCFWorkflowFiltering:
    def __init__(
        self,
        bam_folder,
        ref_genome,
        gene_name,
        gene_positions,
        minipileup_path,
        output_folder,
        threads,
        max_reads
    ):

        self.bam_folder = bam_folder
        self.ref_genome = ref_genome
        self.gene_name = gene_name
        self.gene_positions = gene_positions  # Liste des positions et variantes d'intérêt
        self.minipileup_path = minipileup_path
        self.output_folder = output_folder
        self.threads = threads
        self.max_reads = max_reads  # <--- Paramètre max_reads par défaut 5000

        # Créer les dossiers nécessaires
        self.create_folder(self.output_folder)
        self.per_read_bam_folder = os.path.join(self.output_folder, "per_read_bams")
        self.per_read_vcf_folder = os.path.join(self.output_folder, "per_read_vcfs")
        self.normalized_vcf_folder = os.path.join(self.output_folder, "normalized_vcfs")
        self.create_folder(self.per_read_bam_folder)
        self.create_folder(self.per_read_vcf_folder)
        self.create_folder(self.normalized_vcf_folder)

        # Obtenir les noms des variantes et la liste des barcodes
        self.variant_names = self.get_variant_names()
        self.barcode_list = self.get_barcode_list()

        # Préparer les structures de données
        self.read_variant_records = []  # Stocke les variantes détectées par lecture
        self.total_reads_per_barcode = {}  # Sera rempli avec calculate_total_reads_per_barcode()

        # Fichier read_variants_table.csv
        self.read_variants_path = os.path.join(self.output_folder, 'read_variants_table.csv')
        if not os.path.exists(self.read_variants_path):
            df_init = pd.DataFrame(columns=["Read", "Variants"])
            df_init.to_csv(self.read_variants_path, index=False)
            print(f"Initialisé {self.read_variants_path} avec les en-têtes.")

        # Fichier filtered_reads_count.csv pour enregistrer le nombre de reads uniques filtrés
        self.filtered_reads_count_path = os.path.join(self.output_folder, 'filtered_reads_count.csv')

    def create_folder(self, folder_path):
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            print(f"Créé le dossier : {folder_path}")

    def get_variant_names(self):
        all_variant_names = {
            # CYTB
            ('Cytb', 919, 'C', 'T'): 'S34L',
            ('Cytb', 1205, 'T', 'A'): 'F129L_A',
            ('Cytb', 1205, 'T', 'G'): 'F129L_G',
            ('Cytb', 1246, 'G', 'C'): 'G143A',
            ('Cytb', 1420, 'T', 'C'): 'L201S',
            ('Cytb', 1422, 'C', 'CATGAGG'): 'E203-DE-V204',
            ('Cytb', 1424, 'T', 'TGAGGTA'): 'E203-VE-V204',

            # OSBP
            ('OSBP', 4169, 'T', 'G'): 'L863W',
            ('OSBP', 4091, 'A', 'T'): 'N837I',
            ('OSBP', 3890, 'G', 'T'): 'G770V',

            # TUBULIN
            ('Tubulin', 1541, 'T', 'A'): 'C239S_AGC',
            ('Tubulin', 1541, 'T', 'G'): 'C239G_GGC',
            ('Tubulin', 1542, 'G', 'C'): 'C239S_TCC',

            # ATPase
            ('ATPase', 2857, 'G', 'A'): 'G759S',
        }
        # Filtrer uniquement les variantes pour le gène sélectionné
        return {k: v for k, v in all_variant_names.items() if k[0] == self.gene_name}

    def get_barcode_list(self):
        """Récupère la liste des barcodes à partir des noms de fichiers BAM."""
        barcodes = set()
        for bam_file in os.listdir(self.bam_folder):
            if bam_file.endswith('.bam'):
                barcode = os.path.splitext(bam_file)[0].split('-')[0]
                barcodes.add(barcode)
        return sorted(list(barcodes))


    # F I L T R A G E   ( U N I Q U E S   &   L I M I T É )   D E S   R E A D S

    def filter_unique_reads(self, bam_path):
        unique_reads = set()
        filtered_reads = []
        with pysam.AlignmentFile(bam_path, "rb") as bamfile:
            for read in bamfile.fetch(until_eof=True):
                if read.query_name not in unique_reads:
                    unique_reads.add(read.query_name)
                    filtered_reads.append(read)
                if len(filtered_reads) >= self.max_reads:
                    break
        return filtered_reads, len(unique_reads)

    def save_filtered_reads_count(self, barcode, count):
        df = pd.DataFrame([[barcode, count]], columns=["Barcode", "Filtered_Read_Count"])
        if os.path.exists(self.filtered_reads_count_path):
            df.to_csv(self.filtered_reads_count_path, mode='a', header=False, index=False)
        else:
            df.to_csv(self.filtered_reads_count_path, index=False)
        print(f"Sauvegardé le nombre de reads uniques filtrés pour {barcode}: {count}")


    # B A M   →   P A R   L E C T U R E

    def split_bam_by_read(self):
        """
        Sépare chaque fichier BAM en plusieurs fichiers BAM,
        chacun contenant une seule lecture filtrée.
        """
        bam_files = [f for f in os.listdir(self.bam_folder) if f.endswith('.bam')]
        bam_paths = [os.path.join(self.bam_folder, bam_file) for bam_file in bam_files]
        with Pool(processes=self.threads) as pool:
            pool.map(self.process_bam_file, bam_paths)

    def process_bam_file(self, bam_path):
        """
        Divise un seul fichier BAM en fichiers BAM par lecture,
        en supprimant les doublons et en limitant à self.max_reads.
        """
        bam_file = os.path.basename(bam_path)
        bam_basename = os.path.splitext(bam_file)[0]
        barcode = bam_basename.split('-')[0]
        barcode_folder = os.path.join(self.per_read_bam_folder, barcode)
        self.create_folder(barcode_folder)
        print(f"Traitement du fichier BAM : {bam_file}")

        try:
            filtered_reads, unique_count = self.filter_unique_reads(bam_path)

            # Écrire un sous-fichier BAM par read filtré
            for idx, read in enumerate(filtered_reads, start=1):
                per_read_bam_path = os.path.join(barcode_folder, f"{bam_basename}_read{idx}.bam")
                with pysam.AlignmentFile(per_read_bam_path, "wb", template=pysam.AlignmentFile(bam_path, "rb")) as out_bam:
                    out_bam.write(read)
                pysam.index(per_read_bam_path)

            print(f"Divisé {bam_file} en {len(filtered_reads)} BAMs uniques/filtrés.")
            self.save_filtered_reads_count(barcode, unique_count)
        except Exception as e:
            print(f"Erreur lors du traitement de {bam_file} : {e}")


    #  A P P E L   D E S   V A R I A N T S   P A R   L E C T U R E

    def call_variants_per_read(self):
        bam_folders = [os.path.join(self.per_read_bam_folder, barcode) for barcode in self.barcode_list]
        with Pool(processes=self.threads) as pool:
            pool.map(self.process_bam_folder, bam_folders)
        print("Appel des variantes sur les fichiers BAM par lecture complété.")

    def process_bam_folder(self, barcode_bam_folder):
        barcode = os.path.basename(barcode_bam_folder)
        barcode_vcf_folder = os.path.join(self.per_read_vcf_folder, barcode)
        self.create_folder(barcode_vcf_folder)
        bam_files = [f for f in os.listdir(barcode_bam_folder) if f.endswith('.bam')]

        for bam_file in bam_files:
            bam_path = os.path.join(barcode_bam_folder, bam_file)
            vcf_file = bam_file.replace('.bam', '.vcf')
            vcf_path = os.path.join(barcode_vcf_folder, vcf_file)

            cmd = f"{self.minipileup_path} -vc -f {self.ref_genome} {bam_path}"
            try:
                with open(vcf_path, 'w') as vcf_output:
                    subprocess.run(cmd, shell=True, check=True, stdout=vcf_output, stderr=subprocess.PIPE)
                print(f"Appel des variantes pour {bam_file}, sortie : {vcf_file}")
            except subprocess.CalledProcessError as e:
                print(f"Minipileup a échoué pour {bam_file} : {e.stderr.decode().strip()}")


    # N O R M A L I S A T I O N   D E S   V C F   ( B C F T O O L S )

    def normalize_vcfs(self):
        """Normalise les fichiers VCF via bcftools norm."""
        vcf_folders = [os.path.join(self.per_read_vcf_folder, barcode) for barcode in self.barcode_list]
        with Pool(processes=self.threads) as pool:
            pool.map(self.process_vcf_folder, vcf_folders)
        print("Normalisation des fichiers VCF par lecture complétée.")

    def process_vcf_folder(self, barcode_vcf_folder):
        """Normalise tous les fichiers VCF dans un dossier de barcode via bcftools norm."""
        barcode = os.path.basename(barcode_vcf_folder)
        norm_folder = os.path.join(self.normalized_vcf_folder, barcode)
        self.create_folder(norm_folder)

        vcf_files = [f for f in os.listdir(barcode_vcf_folder) if f.endswith('.vcf')]
        for vcf_file in vcf_files:
            input_vcf = os.path.join(barcode_vcf_folder, vcf_file)
            step1_vcf = os.path.join(norm_folder, f"{os.path.splitext(vcf_file)[0]}_norm_step1.vcf")
            step2_vcf = os.path.join(norm_folder, f"{os.path.splitext(vcf_file)[0]}_norm_step2.vcf")

            try:
                subprocess.run(['bcftools', 'norm', '-m-', '-o', step1_vcf, input_vcf], check=True)
                subprocess.run(['bcftools', 'norm', '-f', self.ref_genome, '-cw', '-o', step2_vcf, step1_vcf], check=True)
                print(f"VCF normalisé : {step2_vcf}")
            except subprocess.CalledProcessError as e:
                print(f"bcftools norm a échoué pour {vcf_file} (barcode {barcode}) : {e.stderr.decode().strip()}")


    # C O L L E C T E   D E S   V A R I A N T S   &   G É N É R A T I O N

    def extract_variants_from_vcf(self, vcf_path, variant_names):

        variants_present = []

        # Variables pour accumuler le alt_count sur les positions "spéciales" E203-DE-V204
        alt_count_1422 = 0
        alt_count_1428 = 0
        alt_count_1422_bis = 0
        alt_count_1424 = 0

        try:
            with open(vcf_path, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    cols = line.strip().split('\t')
                    chrom = cols[0]
                    if chrom != self.gene_name:
                        continue

                    pos = int(cols[1])
                    ref = cols[3]
                    alt = cols[4]

                    # 1) Vérification standard : (chrom, pos, ref, alt) est-il dans variant_names ?
                    key = (chrom, pos, ref, alt)
                    if key in variant_names:
                        variants_present.append(variant_names[key])

                    # 2) Gestion spéciale E203-DE-V204 sur Cytb
                    if self.gene_name == 'Cytb':
                        fmt = cols[8].split(':')
                        sample_vals = cols[9].split(':')
                        if 'AD' in fmt:
                            ad_index = fmt.index('AD')
                            ad_value = sample_vals[ad_index]
                            parts = ad_value.split(',')
                            if len(parts) >= 2 and all(x.isdigit() for x in parts[:2]):
                                alt_count_here = int(parts[1])  # On n'accumule que le alt_count

                                # On ne l'incrémente QUE si (pos, ref, alt) correspond
                                if pos == 1422 and ref == 'C' and alt == 'CATGAGG':
                                    alt_count_1422 += alt_count_here
                                elif pos == 1428 and ref == 'G' and alt == 'GA':
                                    alt_count_1428 += alt_count_here
                                elif pos == 1422 and ref == 'C' and alt == 'CATGAAG':
                                    alt_count_1422_bis += alt_count_here
                                elif pos == 1424 and ref == 'T' and alt == 'TGAGGAC':
                                    alt_count_1424 += alt_count_here
            # À la fin du fichier, si au moins une de ces positions a alt_count > 0,
            # on déclare la variante combinée "E203-DE-V204".
            if self.gene_name == 'Cytb':
                if (alt_count_1422 > 0
                    or alt_count_1428 > 0
                    or alt_count_1422_bis > 0
                    or alt_count_1424 > 0):
                    if "E203-DE-V204" not in variants_present:
                        variants_present.append("E203-DE-V204")

        except Exception as e:
            print(f"Erreur lors de la lecture du fichier {vcf_path} : {e}")

        return variants_present

    def collect_variants_and_generate_files(self):

        for barcode in self.barcode_list:
            barcode_norm_vcf_folder = os.path.join(self.normalized_vcf_folder, barcode)
            if not os.path.exists(barcode_norm_vcf_folder):
                # Il se peut qu'il n'y ait pas eu de BAM (ou autre problème)
                continue

            vcf_files = [f for f in os.listdir(barcode_norm_vcf_folder) if f.endswith('_norm_step2.vcf')]
            total_processed = 0

            for vcf_file in vcf_files:
                vcf_path = os.path.join(barcode_norm_vcf_folder, vcf_file)
                read_variants = self.extract_variants_from_vcf(vcf_path, self.variant_names)
                total_processed += 1

                read_basename = os.path.splitext(vcf_file)[0]
                read_id = read_basename.replace('_norm_step2', '')
                self.read_variant_records.append({'Read': read_id, 'Variants': read_variants})

            self.total_reads_per_barcode[barcode] = total_processed
            print(f"Variantes collectées pour {barcode}, total lectures : {total_processed}")

        # Sauvegarder le tableau read-variants
        self.save_read_variants_table()

        # Générer known_variants_filtered
        self.generate_known_variants_filtered()

        # Générer unknown_variants
        self.generate_unknown_variants()

    def save_read_variants_table(self):
        """Sauvegarde le tableau des lectures et leurs variantes dans un CSV."""
        df = pd.DataFrame(self.read_variant_records)
        if not df.empty:
            df['Variants'] = df['Variants'].apply(lambda vs: str(vs) if vs else '[]')
            df.to_csv(self.read_variants_path, index=False)
            print(f"Tableau Read-Variants sauvegardé : {self.read_variants_path}")
        else:
            print("Aucune donnée de variantes à sauvegarder dans read_variants_table.csv")


    # A U T R E S   F O N C T I O N S

    def calculate_total_reads_per_barcode(self):
        """Calcule (ou recalcule) le nombre total de lectures pour chaque barcode depuis read_variants_table.csv."""
        total_reads = {}
        if os.path.exists(self.read_variants_path):
            df = pd.read_csv(self.read_variants_path)
            if not df.empty and 'Read' in df.columns:
                df['Barcode'] = df['Read'].apply(lambda x: x.split('_')[0].split('-')[0])
                total_reads = df.groupby('Barcode').size().to_dict()
        return total_reads

    def generate_known_variants_filtered(self):

        print("Génération de known_variants_filtered.tsv à partir de read_variants_table.csv...")

        if not os.path.exists(self.read_variants_path):
            print(f"read_variants_table.csv non trouvé : {self.read_variants_path}")
            return

        read_variants_df = pd.read_csv(self.read_variants_path)
        if read_variants_df.empty:
            print("Aucune donnée dans read_variants_table.csv. known_variants_filtered.tsv non généré.")
            return

        if 'Read' not in read_variants_df.columns or 'Variants' not in read_variants_df.columns:
            print("Colonnes attendues manquantes dans read_variants_table.csv.")
            return

        read_variants_df['Barcode'] = read_variants_df['Read'].apply(lambda x: x.split('_')[0].split('-')[0])
        read_variants_df['Variants'] = read_variants_df['Variants'].apply(
            lambda x: ast.literal_eval(x) if isinstance(x, str) else []
        )

        # On déclare comme "connus" : ceux de self.variant_names + "E203-DE-V204"
        known_variants = set(self.variant_names.values()) | {"E203-DE-V204"}

        # On ne garde que les variantes dans known_variants
        read_variants_df['Known_Variants'] = read_variants_df['Variants'].apply(
            lambda vs: [v for v in vs if v in known_variants]
        )

        known_exploded = read_variants_df.explode('Known_Variants')
        known_exploded = known_exploded.dropna(subset=['Known_Variants'])

        # Count occurrences par (Barcode, Known_Variants)
        known_counts = known_exploded.groupby(['Barcode', 'Known_Variants']).size().reset_index(name='Count')
        pivot_counts = known_counts.pivot_table(
            index='Known_Variants',
            columns='Barcode',
            values='Count',
            fill_value=0
        ).reset_index()
        pivot_counts.rename(columns={'Known_Variants': 'variant_name'}, inplace=True)

        # Recalcule total_reads_per_barcode
        self.total_reads_per_barcode = self.calculate_total_reads_per_barcode()

        # Ajouter {barcode}_VARIANT_COUNT et {barcode}_VARIANT_FREQ
        for barcode in self.barcode_list:
            count_col = f"{barcode}_VARIANT_COUNT"
            freq_col = f"{barcode}_VARIANT_FREQ"

            total_reads = self.total_reads_per_barcode.get(barcode, 1)
            if barcode in pivot_counts.columns:
                pivot_counts[freq_col] = pivot_counts[barcode].apply(
                    lambda c: f"{(c / total_reads * 100):.2f}%" if total_reads > 0 else "0.00%"
                )
                pivot_counts.rename(columns={barcode: count_col}, inplace=True)
            else:
                pivot_counts[count_col] = 0
                pivot_counts[freq_col] = "0.00%"

        cols_order = ['variant_name']
        for bc in self.barcode_list:
            cols_order += [f"{bc}_VARIANT_COUNT", f"{bc}_VARIANT_FREQ"]
        known_variants_filtered_df = pivot_counts[cols_order]

        out_path = os.path.join(self.output_folder, 'known_variants_filtered.tsv')
        known_variants_filtered_df.to_csv(out_path, sep='\t', index=False)
        print(f"known_variants_filtered.tsv sauvegardé : {out_path}")

    def generate_unknown_variants(self):
        print("Génération de unknown_variants.tsv via minipileup sur fichiers BAM entiers...")

        unknown_variants_dict = {}

        known_path = os.path.join(self.output_folder, 'known_variants_filtered.tsv')
        if not os.path.exists(known_path):
            print(f"Fichier known_variants_filtered.tsv introuvable : {known_path}")
            return

        known_variants_df = pd.read_csv(known_path, sep='\t')
        if 'variant_name' not in known_variants_df.columns:
            print("Colonne 'variant_name' introuvable dans known_variants_filtered.tsv.")
            return

        # On récupère la liste de toutes les "variante_name" déjà cataloguées
        known_variants_list = set(known_variants_df['variant_name'])

        bam_files = [f for f in os.listdir(self.bam_folder) if f.endswith('.bam')]
        for bam_file in bam_files:
            bam_path = os.path.join(self.bam_folder, bam_file)
            print(f"Minipileup sur : {bam_file}")
            try:
                temp_vcf = os.path.join(self.output_folder, f"{bam_file}_temp.vcf")
                cmd = f"{self.minipileup_path} -vc -f {self.ref_genome} {bam_path}"
                with open(temp_vcf, 'w') as out:
                    subprocess.run(cmd, shell=True, check=True, stdout=out, stderr=subprocess.PIPE)

                # Normalisation
                norm_temp_vcf = os.path.join(self.output_folder, f"{bam_file}_temp_norm.vcf")
                try:
                    subprocess.run(['bcftools', 'norm', '-m-', '-o', norm_temp_vcf, temp_vcf], check=True)
                    subprocess.run(['bcftools', 'norm', '-f', self.ref_genome, '-cw', '-o', temp_vcf, norm_temp_vcf], check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Erreur bcftools norm sur {temp_vcf} : {e.stderr.decode().strip()}")
                    continue

                # Parse ce VCF
                variants = self.parse_minipileup_output_file(temp_vcf)
                barcode = bam_file.split('-')[0]

                for variant_info in variants:
                    variant = variant_info['variant']
                    ad = variant_info['AD']

                    # On parse 'CHROM:POS:REF>ALT'
                    try:
                        chrom, pos_str, ref_alt = variant.split(':')
                        pos = int(pos_str)
                        ref, alt = ref_alt.split('>')
                    except Exception as e:
                        print(f"Erreur parsing variant '{variant}': {e}")
                        continue

                    # On tente de voir si c'est connu via get_variant_name
                    variant_name = self.get_variant_name(chrom, pos, ref, alt)
                    if variant_name in known_variants_list:
                        continue  # c'est connu, on ignore

                    # Ajout dans unknown_variants_dict
                    key = (chrom, pos, ref, alt)
                    if key not in unknown_variants_dict:
                        unknown_variants_dict[key] = {'CHROM': chrom, 'POS': pos, 'REF': ref, 'ALT': alt}

                    try:
                        ref_count, alt_count = map(int, ad.split(','))
                    except:
                        ref_count, alt_count = 0, 0

                    total_count = ref_count + alt_count
                    frequency = (alt_count / total_count * 100) if total_count > 0 else 0.0
                    unknown_variants_dict[key][barcode] = {
                        'AD': ad,
                        'Frequency': f"{frequency:.2f}%",
                        'Count': alt_count
                    }

                os.remove(temp_vcf)
                os.remove(norm_temp_vcf)
            except subprocess.CalledProcessError as e:
                print(f"Minipileup échoué sur {bam_file} : {e.stderr.decode().strip()}")

        if not unknown_variants_dict:
            print("Aucune variante inconnue détectée.")
            return

        # Construction DataFrame final
        unknown_variants_list = []
        barcodes_sorted = self.barcode_list
        for key, data in unknown_variants_dict.items():
            row = [data['CHROM'], data['POS'], data['REF'], data['ALT']]
            for bc in barcodes_sorted:
                if bc in data:
                    row += [data[bc]['AD'], data[bc]['Frequency'], data[bc]['Count']]
                else:
                    row += ['0,0', '0.00%', '0']
            unknown_variants_list.append(row)

        header = ['CHROM', 'POS', 'REF', 'ALT']
        for bc in barcodes_sorted:
            header += [f"{bc}_AD", f"{bc}_VARIANT_FREQ", f"{bc}_VARIANT_COUNT"]
        df_unknown = pd.DataFrame(unknown_variants_list, columns=header)

        out_unknown_path = os.path.join(self.output_folder, 'unknown_variants.tsv')
        df_unknown.to_csv(out_unknown_path, sep='\t', index=False)
        print(f"unknown_variants.tsv sauvegardé : {out_unknown_path}")

    def get_variant_name_from_variant_string(self, variant_string):
        try:
            chrom, pos, ref_alt = variant_string.split(':')
            ref, alt = ref_alt.split('>')
            return self.get_variant_name(chrom, pos, ref, alt)
        except:
            return None

    def get_variant_name(self, chrom, pos, ref, alt):
        try:
            pos_int = int(pos)
        except:
            return None
        key = (chrom, pos_int, ref, alt)
        return self.variant_names.get(key, None)

    def parse_minipileup_output_file(self, vcf_path):
        variants = []
        try:
            with open(vcf_path, 'r') as vcf:
                for line in vcf:
                    if line.startswith('#'):
                        continue
                    cols = line.strip().split('\t')
                    chrom = cols[0]
                    pos = int(cols[1])
                    ref = cols[3]
                    alt = cols[4]
                    fmt = cols[8].split(':')
                    sample_vals = cols[9].split(':')

                    if 'AD' in fmt:
                        ad_index = fmt.index('AD')
                        ad_value = sample_vals[ad_index]
                        parts = ad_value.split(',')
                        if len(parts) >= 2 and all(x.isdigit() for x in parts[:2]):
                            ref_count = int(parts[0])
                            alt_count = int(parts[1])
                            ad_value = f"{ref_count},{alt_count}"
                        else:
                            ad_value = "0,0"
                    else:
                        ad_value = "0,0"

                    variant_id = f"{chrom}:{pos}:{ref}>{alt}"
                    variants.append({'variant': variant_id, 'AD': ad_value})
        except Exception as e:
            print(f"Erreur lecture {vcf_path} : {e}")
        return variants


    # C A L C U L   D E S   F R É Q U E N C E S   W I L D - T Y P E

    def calculate_wildtype_frequencies(self):
        out_file = os.path.join(self.output_folder, "wildtype_frequencies.csv")
        results = []
        for barcode in self.barcode_list:
            total_reads = self.total_reads_per_barcode.get(barcode, 0)
            wildtype_count = 0
            for rec in self.read_variant_records:
                if rec['Read'].startswith(barcode):
                    if not rec['Variants']:
                        wildtype_count += 1
            freq_wt = (wildtype_count / total_reads * 100) if total_reads > 0 else 0.0
            freq_wt = round(freq_wt, 2)
            results.append([barcode, wildtype_count, f"{freq_wt}%"])
            print(f"Barcode : {barcode} | Total : {total_reads} | Wild-type : {wildtype_count} | Freq : {freq_wt}%")

        df_wt = pd.DataFrame(results, columns=['Sample', 'Count', 'Frequency'])
        # Tri par numéro de barcode (optionnel)
        def extract_barcode_num(x):
            m = re.search(r'barcode(\d+)', x)
            return int(m.group(1)) if m else 999999
        df_wt['Barcode_Num'] = df_wt['Sample'].apply(extract_barcode_num)
        df_wt = df_wt.sort_values('Barcode_Num').drop('Barcode_Num', axis=1)

        df_wt.to_csv(out_file, index=False)
        print(f"Fréquences wild-type sauvegardées : {out_file}")

    #
    # N E T T O Y A G E   D E S   F I C H I E R S   T E M P
    #
    def clean_up(self):
        """Supprime les dossiers intermédiaires (BAM par lecture, VCF par lecture, VCF normalisés)."""
        dirs = [self.per_read_bam_folder, self.per_read_vcf_folder, self.normalized_vcf_folder]
        for d in dirs:
            if os.path.exists(d):
                shutil.rmtree(d)
                print(f"Supprimé : {d}")
            else:
                print(f"Dossier introuvable (ignoré) : {d}")



    def run(self):
        """Exécute le workflow complet."""
        # Étape 1 : BAM → par lecture
        self.split_bam_by_read()

        # Étape 2 : Appel variants minipileup par lecture
        self.call_variants_per_read()

        # Étape 3 : Normalisation bcftools
        self.normalize_vcfs()

        # Étape 4 : Collecte + Génération read_variants_table / known / unknown
        self.collect_variants_and_generate_files()

        # Étape 5 : Fréquences wild-type
        self.calculate_wildtype_frequencies()

        # Étape 6 : Nettoyage
        self.clean_up()



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

bam_folder = sys.argv[1]
ref_genome = sys.argv[2]
gene_name = str(sys.argv[3])
minipileup_path = sys.argv[4]
output_folder = sys.argv[5]
threads = int(sys.argv[6])
max_reads = int(sys.argv[7])

workflow = VCFWorkflowFiltering(
    bam_folder=bam_folder,
    ref_genome=ref_genome,
    gene_name=gene_name,
    gene_positions=gene_positions.get(gene_name, []),
    minipileup_path=minipileup_path,
    output_folder=output_folder,
    threads=threads,
    max_reads=max_reads
)
workflow.run()
