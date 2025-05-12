import sys
import os
import pandas as pd
import subprocess
import numpy as np
import re
import pysam
import shutil
from multiprocessing import Pool
import ast

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
        self.gene_positions = gene_positions
        self.minipileup_path = minipileup_path
        self.output_folder = output_folder
        self.threads = threads
        self.max_reads = max_reads

        # Créer les dossiers nécessaires
        self.create_folder(self.output_folder)
        self.per_read_bam_folder = os.path.join(self.output_folder, "per_read_bams")
        self.per_read_vcf_folder = os.path.join(self.output_folder, "per_read_vcfs")
        self.normalized_vcf_folder = os.path.join(self.output_folder, "normalized_vcfs")
        self.create_folder(self.per_read_bam_folder)
        self.create_folder(self.per_read_vcf_folder)
        self.create_folder(self.normalized_vcf_folder)

        # Obtenir les noms des variantes (nouveau dictionnaire)
        self.variant_names = self.get_variant_names()
        self.barcode_list = self.get_barcode_list()

        # Données
        self.read_variant_records = []
        self.total_reads_per_barcode = {}

        # read_variants_table.csv
        self.read_variants_path = os.path.join(
            self.output_folder, "read_variants_table.csv"
        )
        if not os.path.exists(self.read_variants_path):
            df_init = pd.DataFrame(columns=["Read", "Variants"])
            df_init.to_csv(self.read_variants_path, index=False)
            print(f"Initialisé {self.read_variants_path} avec les en-têtes.")

        # filtered_reads_count.csv
        self.filtered_reads_count_path = os.path.join(
            self.output_folder, "filtered_reads_count.csv"
        )

    def create_folder(self, folder_path):
        """Crée un dossier s'il n'existe pas déjà."""
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            print(f"Créé le dossier : {folder_path}")

    # AJOUT C136

    @staticmethod
    def consolidate_c136(variants_list):
        """
        Remplace C136_uno + C136_dos par C136 quand ils co-existent sur un même read.
        """
        if "C136_uno" in variants_list and "C136_dos" in variants_list:
            return [v for v in variants_list if v not in ("C136_uno", "C136_dos")] + [
                "C136"
            ]
        return [v for v in variants_list if v not in ("C136_uno", "C136_dos")]

    # ------------------------------------------------------------------

    def get_variant_names(self):
        all_variant_names = {
            # CYP51
            ("CYP51", 149, "T", "C"): "L50S",
            ("CYP51", 320, "A", "T"): "D107V",
            ("CYP51", 401, "A", "G"): "D134G",
            ("CYP51", 407, "T", "C"): "V136A",
            ("CYP51", 410, "A", "T"): "Y137F",
            ("CYP51", 433, "A", "C"): "M145L_CTG",
            ("CYP51", 433, "A", "T"): "M145L_TTG",
            ("CYP51", 590, "A", "G"): "N178S",
            ("CYP51", 620, "G", "A"): "S188N",
            ("CYP51", 679, "T", "A"): "S208T",
            ("CYP51", 841, "T", "C"): "Y262H",
            ("CYP51", 1398, "C", "T"): "R376C",
            ("CYP51", 1408, "C", "G"): "A379G",
            ("CYP51", 1413, "A", "G"): "I381V",
            ("CYP51", 1648, "A", "C"): "Y459S",
            ("CYP51", 1645, "ACTATGG", "A"): "ΔY459/G460",
            ("CYP51", 1654, "A", "C"): "Y461S",
            ("CYP51", 1653, "T", "C"): "Y461H",
            ("CYP51", 1811, "C", "A"): "N513K_AAA",
            ("CYP51", 1811, "C", "G"): "N513K_AAG",
            ("CYP51", 1843, "G", "C"): "S524T",
            ("CYP51", 1500, "G", "A"): "A410T",
            # --- Composite C136 (uno / dos)
            ("CYP51", 406, "G", "T"): "C136_uno",
            ("CYP51", 407, "T", "G"): "C136_dos",
            # CytB
            ("CytB", 110, "G", "T"): "G37V",
            ("CytB", 428, "G", "C"): "G143A",
            # SdhB
            ("SdhB", 788, "A", "C"): "N225T",
            ("SdhB", 788, "A", "T"): "N225I",
            ("SdhB", 908, "G", "C"): "R265P",
            ("SdhB", 914, "A", "T"): "H267L",
            ("SdhB", 917, "C", "T"): "T268I",
            ("SdhB", 916, "A", "G"): "T268A",
            # SdhC
            ("SdhC", 365, "C", "A"): "T79N",
            ("SdhC", 368, "G", "C"): "W80S",
            ("SdhC", 379, "G", "T"): "A84F",
            ("SdhC", 386, "A", "G"): "N86S",
            ("SdhC", 387, "C", "A"): "N86K_AAA",
            ("SdhC", 387, "C", "G"): "N86K_AAG",
            ("SdhC", 508, "C", "G"): "P127A",
            ("SdhC", 582, "G", "C"): "R151S_AGC",
            ("SdhC", 582, "G", "T"): "R151S_AGT",
            ("SdhC", 581, "G", "C"): "R151T",
            ("SdhC", 581, "G", "T"): "R151M",
            ("SdhC", 580, "A", "G"): "R151G",
            ("SdhC", 584, "A", "G"): "H152R",
            ("SdhC", 611, "T", "G"): "I161S",
            ("SdhC", 625, "G", "A"): "V166M",
            ("SdhC", 632, "C", "G"): "T168R",
            # SdhD
            ("SdhD", 139, "C", "T"): "R47W",
            ("SdhD", 148, "A", "T"): "I50F",
            ("SdhD", 148, "A", "C"): "I50L",
            ("SdhD", 399, "A", "G"): "M114V",
            ("SdhD", 497, "A", "G"): "D129G",
            ("SdhD", 498, "C", "A"): "D129E_GAA",
            ("SdhD", 498, "C", "G"): "D129E_GAG",
        }
        # Filtrer uniquement le gène choisi
        return {k: v for k, v in all_variant_names.items() if k[0] == self.gene_name}

    def get_barcode_list(self):
        """Récupère la liste des barcodes à partir des noms de fichiers BAM."""
        barcodes = set()
        for bam_file in os.listdir(self.bam_folder):
            if bam_file.endswith(".bam"):
                barcode = os.path.splitext(bam_file)[0].split("-")[0]
                barcodes.add(barcode)
        return sorted(list(barcodes))

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
            df.to_csv(
                self.filtered_reads_count_path, mode="a", header=False, index=False
            )
        else:
            df.to_csv(self.filtered_reads_count_path, index=False)
        print(f"Sauvegardé le nombre de reads uniques filtrés pour {barcode}: {count}")

    def split_bam_by_read(self):
        bam_files = [f for f in os.listdir(self.bam_folder) if f.endswith(".bam")]
        bam_paths = [os.path.join(self.bam_folder, bf) for bf in bam_files]
        with Pool(processes=self.threads) as pool:
            pool.map(self.process_bam_file, bam_paths)

    def process_bam_file(self, bam_path):
        bam_file = os.path.basename(bam_path)
        bam_basename = os.path.splitext(bam_file)[0]
        barcode = bam_basename.split("-")[0]
        barcode_folder = os.path.join(self.per_read_bam_folder, barcode)
        self.create_folder(barcode_folder)

        print(f"Traitement du BAM: {bam_file}")
        try:
            filtered_reads, unique_count = self.filter_unique_reads(bam_path)
            for idx, read in enumerate(filtered_reads, start=1):
                out_bam_path = os.path.join(
                    barcode_folder, f"{bam_basename}_read{idx}.bam"
                )
                with pysam.AlignmentFile(
                    out_bam_path, "wb", template=pysam.AlignmentFile(bam_path, "rb")
                ) as out_bam:
                    out_bam.write(read)
                pysam.index(out_bam_path)
            print(f"Divisé {bam_file} => {len(filtered_reads)} BAM uniques/filtrés.")
            self.save_filtered_reads_count(barcode, unique_count)
        except Exception as e:
            print(f"Erreur sur {bam_file}: {e}")

    # ------------------------------------------------------------------
    # MINIPILEUP PAR LECTURE (inchangé)
    # ------------------------------------------------------------------
    def call_variants_per_read(self):
        bam_folders = [
            os.path.join(self.per_read_bam_folder, bc) for bc in self.barcode_list
        ]
        with Pool(processes=self.threads) as pool:
            pool.map(self.process_bam_folder, bam_folders)
        print("Appel variants minipileup par lecture terminé.")

    def process_bam_folder(self, barcode_bam_folder):
        barcode = os.path.basename(barcode_bam_folder)
        barcode_vcf_folder = os.path.join(self.per_read_vcf_folder, barcode)
        self.create_folder(barcode_vcf_folder)

        bam_files = [
            f for f in os.listdir(barcode_bam_folder) if f.endswith(".bam")
        ]
        for bam_file in bam_files:
            bam_path = os.path.join(barcode_bam_folder, bam_file)
            vcf_file = bam_file.replace(".bam", ".vcf")
            vcf_path = os.path.join(barcode_vcf_folder, vcf_file)

            cmd = f"{self.minipileup_path} -vc -f {self.ref_genome} {bam_path}"
            try:
                with open(vcf_path, "w") as out_vcf:
                    subprocess.run(
                        cmd, shell=True, check=True, stdout=out_vcf, stderr=subprocess.PIPE
                    )
                print(f"Minipileup => {bam_file}, output: {vcf_file}")
            except subprocess.CalledProcessError as e:
                print(f"Minipileup échoué: {e.stderr.decode().strip()}")

    # ------------------------------------------------------------------
    # NORMALISATION DES VCF (inchangé)
    # ------------------------------------------------------------------
    def normalize_vcfs(self):
        vcf_folders = [
            os.path.join(self.per_read_vcf_folder, bc) for bc in self.barcode_list
        ]
        with Pool(processes=self.threads) as pool:
            pool.map(self.process_vcf_folder, vcf_folders)
        print("Normalisation des VCF complétée.")

    def process_vcf_folder(self, barcode_vcf_folder):
        barcode = os.path.basename(barcode_vcf_folder)
        norm_folder = os.path.join(self.normalized_vcf_folder, barcode)
        self.create_folder(norm_folder)

        vcf_files = [f for f in os.listdir(barcode_vcf_folder) if f.endswith(".vcf")]
        for vcf_file in vcf_files:
            input_vcf = os.path.join(barcode_vcf_folder, vcf_file)
            step1_vcf = os.path.join(
                norm_folder, f"{os.path.splitext(vcf_file)[0]}_norm_step1.vcf"
            )
            step2_vcf = os.path.join(
                norm_folder, f"{os.path.splitext(vcf_file)[0]}_norm_step2.vcf"
            )

            try:
                subprocess.run(["bcftools", "norm", "-m-", "-o", step1_vcf, input_vcf], check=True)
                subprocess.run(
                    ["bcftools", "norm", "-f", self.ref_genome, "-cw", "-o", step2_vcf, step1_vcf],
                    check=True,
                )
                print(f"VCF normalisé => {step2_vcf}")
            except subprocess.CalledProcessError as e:
                print(f"bcftools norm échoué pour {vcf_file} : {e.stderr.decode().strip()}")

    # ------------------------------------------------------------------
    # COLLECTE DES VARIANTS & GÉNÉRATION
    # ------------------------------------------------------------------
    def extract_variants_from_vcf(self, vcf_path, variant_map):
        """
        Lit un VCF normalisé. Pour chaque ligne où CHROM == self.gene_name,
        si (chrom, pos, ref, alt) figure dans variant_map => on l'ajoute.
        """
        variants_present = []
        try:
            with open(vcf_path, "r") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    cols = line.strip().split("\t")
                    if len(cols) < 5:
                        continue
                    chrom = cols[0]
                    if chrom != self.gene_name:
                        continue
                    pos = int(cols[1])
                    ref = cols[3]
                    alt = cols[4]
                    key = (chrom, pos, ref, alt)
                    if key in variant_map:
                        variants_present.append(variant_map[key])
        except Exception as e:
            print(f"Erreur lecture {vcf_path}: {e}")
        return variants_present

    def collect_variants_and_generate_files(self):
        for barcode in self.barcode_list:
            barcode_norm_folder = os.path.join(self.normalized_vcf_folder, barcode)
            if not os.path.exists(barcode_norm_folder):
                continue
            vcf_files = [
                f for f in os.listdir(barcode_norm_folder) if f.endswith("_norm_step2.vcf")
            ]
            total_processed = 0

            for vcf_file in vcf_files:
                vcf_path = os.path.join(barcode_norm_folder, vcf_file)
                read_vars = self.extract_variants_from_vcf(vcf_path, self.variant_names)
                read_vars = self.consolidate_c136(read_vars)  # <-- ajout C136
                total_processed += 1

                read_basename = os.path.splitext(vcf_file)[0]
                read_id = read_basename.replace("_norm_step2", "")
                self.read_variant_records.append({"Read": read_id, "Variants": read_vars})

            self.total_reads_per_barcode[barcode] = total_processed
            print(f"{barcode}: {total_processed} lectures => variants collectés.")

        # Sauvegarde, puis génération des fichiers de sortie
        self.save_read_variants_table()
        self.generate_known_variants_filtered()
        self.generate_unknown_variants()

    def save_read_variants_table(self):
        df = pd.DataFrame(self.read_variant_records)
        if not df.empty:
            df["Variants"] = df["Variants"].apply(lambda vs: str(vs) if vs else "[]")
            df.to_csv(self.read_variants_path, index=False)
            print(f"Tableau read_variants_table => {self.read_variants_path}")
        else:
            print("Aucune variante collectée => read_variants_table non créé.")

    # ------------------------------------------------------------------
    # KNOWN / UNKNOWN / WILDTYPE
    # ------------------------------------------------------------------
    def calculate_total_reads_per_barcode(self):
        totals = {}
        if os.path.exists(self.read_variants_path):
            df = pd.read_csv(self.read_variants_path)
            if not df.empty and "Read" in df.columns:
                df["Barcode"] = df["Read"].apply(
                    lambda x: x.split("_")[0].split("-")[0]
                )
                totals = df.groupby("Barcode").size().to_dict()
        return totals

    def generate_known_variants_filtered(self):
        print("Génération known_variants_filtered.tsv...")
        if not os.path.exists(self.read_variants_path):
            print("read_variants_table.csv introuvable.")
            return
        df = pd.read_csv(self.read_variants_path)
        if df.empty:
            print("read_variants_table.csv vide => known_variants_filtered non généré.")
            return
        if "Read" not in df.columns or "Variants" not in df.columns:
            print("Colonnes manquantes => annulation.")
            return

        df["Barcode"] = df["Read"].apply(lambda x: x.split("_")[0].split("-")[0])
        df["Variants"] = df["Variants"].apply(
            lambda x: ast.literal_eval(x) if isinstance(x, str) else []
        )
        df["Variants"] = df["Variants"].apply(self.consolidate_c136)  # <-- ajout C136

        known_variants = set(self.variant_names.values()) | {"C136"}  # <-- ajout C136
        df["Known_Variants"] = df["Variants"].apply(
            lambda vs: [v for v in vs if v in known_variants]
        )

        known_exploded = df.explode("Known_Variants").dropna(
            subset=["Known_Variants"]
        )
        if known_exploded.empty:
            print("Aucun variant connu => known_variants_filtered vide.")
            out_empty = os.path.join(self.output_folder, "known_variants_filtered.tsv")
            pd.DataFrame(columns=["variant_name"]).to_csv(
                out_empty, sep="\t", index=False
            )
            return

        known_counts = (
            known_exploded.groupby(["Barcode", "Known_Variants"])
            .size()
            .reset_index(name="Count")
        )
        pivot_counts = (
            known_counts.pivot_table(
                index="Known_Variants",
                columns="Barcode",
                values="Count",
                fill_value=0,
            )
            .reset_index()
        )
        pivot_counts.rename(columns={"Known_Variants": "variant_name"}, inplace=True)

        self.total_reads_per_barcode = self.calculate_total_reads_per_barcode()

        for bc in self.barcode_list:
            count_col = f"{bc}_VARIANT_COUNT"
            freq_col = f"{bc}_VARIANT_FREQ"
            tot = self.total_reads_per_barcode.get(bc, 1)
            if bc in pivot_counts.columns:
                pivot_counts[freq_col] = pivot_counts[bc].apply(
                    lambda c: f"{(c / tot * 100):.2f}%" if tot > 0 else "0.00%"
                )
                pivot_counts.rename(columns={bc: count_col}, inplace=True)
            else:
                pivot_counts[count_col] = 0
                pivot_counts[freq_col] = "0.00%"

        cols = ["variant_name"]
        for bc in self.barcode_list:
            cols += [f"{bc}_VARIANT_COUNT", f"{bc}_VARIANT_FREQ"]
        known_variants_filtered_df = pivot_counts[cols]

        out_path = os.path.join(self.output_folder, "known_variants_filtered.tsv")
        known_variants_filtered_df.to_csv(out_path, sep="\t", index=False)
        print(f"known_variants_filtered.tsv => {out_path}")
    def generate_unknown_variants(self):
        print("Génération de unknown_variants.tsv via minipileup sur BAM complet...")
        unknown_dict = {}

        known_path = os.path.join(self.output_folder,'known_variants_filtered.tsv')
        if not os.path.exists(known_path):
            print("known_variants_filtered.tsv introuvable => skip.")
            return

        known_df = pd.read_csv(known_path, sep='\t')
        if 'variant_name' not in known_df.columns:
            print("'variant_name' manquant => skip unknown_variants.")
            return

        known_list = set(known_df['variant_name'])

        bam_files = [f for f in os.listdir(self.bam_folder) if f.endswith('.bam')]
        for bam_file in bam_files:
            print(f"Minipileup => {bam_file}")
            bam_path = os.path.join(self.bam_folder,bam_file)
            try:
                temp_vcf = os.path.join(self.output_folder, f"{bam_file}_temp.vcf")
                cmd = f"{self.minipileup_path} -vc -f {self.ref_genome} {bam_path}"
                with open(temp_vcf, 'w') as out:
                    subprocess.run(cmd, shell=True, check=True, stdout=out, stderr=subprocess.PIPE)

                norm_temp = os.path.join(self.output_folder, f"{bam_file}_temp_norm.vcf")
                try:
                    subprocess.run(['bcftools','norm','-m-','-o', norm_temp, temp_vcf], check=True)
                    subprocess.run(['bcftools','norm','-f', self.ref_genome, '-cw','-o', temp_vcf, norm_temp], check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Erreur bcftools norm: {e.stderr.decode().strip()}")
                    continue

                # Parse
                variants = self.parse_minipileup_output_file(temp_vcf)
                barcode = bam_file.split('-')[0]

                for var_info in variants:
                    variant = var_info['variant']
                    ad = var_info['AD']
                    try:
                        chrom, pos_str, ref_alt = variant.split(':')
                        pos = int(pos_str)
                        ref, alt = ref_alt.split('>')
                    except:
                        continue

                    var_name = self.get_variant_name(chrom,pos,ref,alt)
                    if var_name in known_list:
                        continue  # skip known

                    key = (chrom, pos, ref, alt)
                    if key not in unknown_dict:
                        unknown_dict[key] = {'CHROM': chrom, 'POS':pos, 'REF':ref, 'ALT':alt}
                    try:
                        rc, ac = map(int, ad.split(','))
                    except:
                        rc, ac = 0,0
                    tot = rc+ac
                    freq = (ac/tot*100) if tot>0 else 0
                    unknown_dict[key][barcode] = {
                        'AD': ad,
                        'Frequency': f"{freq:.2f}%",
                        'Count': ac
                    }

                os.remove(temp_vcf)
                os.remove(norm_temp)
            except subprocess.CalledProcessError as e:
                print(f"Minipileup échoué sur {bam_file}: {e.stderr.decode().strip()}")

        if not unknown_dict:
            print("Aucune variante inconnue détectée.")
            return

        # Convert to DataFrame
        unknown_list = []
        for key,data in unknown_dict.items():
            row = [data['CHROM'], data['POS'], data['REF'], data['ALT']]
            for bc in self.barcode_list:
                if bc in data:
                    row += [data[bc]['AD'], data[bc]['Frequency'], data[bc]['Count']]
                else:
                    row += ['0,0','0.00%','0']
            unknown_list.append(row)

        header = ['CHROM','POS','REF','ALT']
        for bc in self.barcode_list:
            header += [f"{bc}_AD", f"{bc}_VARIANT_FREQ", f"{bc}_VARIANT_COUNT"]
        df_unk = pd.DataFrame(unknown_list, columns=header)

        out_unk = os.path.join(self.output_folder,'unknown_variants.tsv')
        df_unk.to_csv(out_unk, sep='\t', index=False)
        print(f"unknown_variants.tsv => {out_unk}")

    def get_variant_name_from_variant_string(self, variant_string):
        """Extrait nom du variant depuis 'CHROM:POS:REF>ALT' si figure dans self.variant_names."""
        try:
            chrom, pos, ref_alt = variant_string.split(':')
            ref, alt = ref_alt.split('>')
            return self.get_variant_name(chrom, pos, ref, alt)
        except:
            return None

    def get_variant_name(self, chrom, pos, ref, alt):
        key = None
        try:
            pos_int = int(pos)
            key = (chrom,pos_int,ref,alt)
        except:
            return None
        return self.variant_names.get(key, None)

    def parse_minipileup_output_file(self, vcf_path):
        """
        Lit un VCF minipileup + AD, renvoie [{'variant':'CHROM:POS:REF>ALT','AD':'x,y'},...].
        """
        variants = []
        try:
            with open(vcf_path,'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    cols = line.strip().split('\t')
                    if len(cols)<10:
                        continue
                    chrom = cols[0]
                    pos = int(cols[1])
                    ref = cols[3]
                    alt = cols[4]
                    fmt = cols[8].split(':')
                    sample_vals = cols[9].split(':')

                    ad_value = "0,0"
                    if 'AD' in fmt:
                        ad_idx = fmt.index('AD')
                        if ad_idx < len(sample_vals):
                            raw_ad = sample_vals[ad_idx]
                            prts = raw_ad.split(',')
                            if len(prts)>=2 and all(x.isdigit() for x in prts[:2]):
                                ad_value = f"{prts[0]},{prts[1]}"

                    variant_id = f"{chrom}:{pos}:{ref}>{alt}"
                    variants.append({'variant': variant_id, 'AD': ad_value})
        except Exception as e:
            print(f"Erreur lecture {vcf_path}: {e}")
        return variants


    # FREQ WILD-TYPE

    def calculate_wildtype_frequencies(self):
        out_file = os.path.join(self.output_folder,"wildtype_frequencies.csv")
        results = []
        for bc in self.barcode_list:
            total_reads = self.total_reads_per_barcode.get(bc,0)
            wt_count = 0
            for rec in self.read_variant_records:
                if rec['Read'].startswith(bc):
                    if not rec['Variants']:
                        wt_count+=1
            freq = round((wt_count/total_reads*100),2) if total_reads>0 else 0.0
            results.append([bc, wt_count, f"{freq}%"])
            print(f"{bc}: total={total_reads}, wildtype={wt_count}, freq={freq}%")

        df_wt = pd.DataFrame(results, columns=['Sample','Count','Frequency'])
        def extract_num(x):
            m = re.search(r'barcode(\d+)', x)
            return int(m.group(1)) if m else 999999
        df_wt['Barcode_Num'] = df_wt['Sample'].apply(extract_num)
        df_wt.sort_values('Barcode_Num', inplace=True)
        df_wt.drop(columns=['Barcode_Num'], inplace=True)

        df_wt.to_csv(out_file, index=False)
        print(f"Fréquences wild-type => {out_file}")


    # NETTOYAGE

    def clean_up(self):
        dirs = [self.per_read_bam_folder,self.per_read_vcf_folder,self.normalized_vcf_folder]
        for d in dirs:
            if os.path.exists(d):
                shutil.rmtree(d)
                print(f"Supprimé: {d}")
            else:
                print(f"Introuvable (ignoré): {d}")


    # RUN

    def run(self):
        """Exécute le workflow complet (read-by-read + minipileup + bcftools)."""
        # 1) BAM -> par lecture
        self.split_bam_by_read()

        # 2) minipileup par lecture
        self.call_variants_per_read()

        # 3) Normalisation
        self.normalize_vcfs()

        # 4) Collect + known/unknown
        self.collect_variants_and_generate_files()

        # 5) Freq wild-type
        self.calculate_wildtype_frequencies()

        # 6) Nettoyage
        self.clean_up()


gene_positions = {
    'CYP51': [
        {'CHROM': 'CYP51', 'POS': 149, 'REF': 'T', 'ALT': 'C'},    # L50S
        {'CHROM': 'CYP51', 'POS': 320, 'REF': 'A', 'ALT': 'T'},    # D107V
        {'CHROM': 'CYP51', 'POS': 401, 'REF': 'A', 'ALT': 'G'},    # D134G
        {'CHROM': 'CYP51', 'POS': 407, 'REF': 'T', 'ALT': 'C'},    # V136A
        {'CHROM': 'CYP51', 'POS': 410, 'REF': 'A', 'ALT': 'T'},    # Y137F
        {'CHROM': 'CYP51', 'POS': 433, 'REF': 'A', 'ALT': 'C'},    # M145L_CTG
        {'CHROM': 'CYP51', 'POS': 433, 'REF': 'A', 'ALT': 'T'},    # M145L_TTG
        {'CHROM': 'CYP51', 'POS': 590, 'REF': 'A', 'ALT': 'G'},    # N178S
        {'CHROM': 'CYP51', 'POS': 620, 'REF': 'G', 'ALT': 'A'},    # S188N
        {'CHROM': 'CYP51', 'POS': 679, 'REF': 'T', 'ALT': 'A'},    # S208T
        {'CHROM': 'CYP51', 'POS': 841, 'REF': 'T', 'ALT': 'C'},    # Y262H
        {'CHROM': 'CYP51', 'POS': 1398, 'REF': 'C', 'ALT': 'T'},   # R376C
        {'CHROM': 'CYP51', 'POS': 1408, 'REF': 'C', 'ALT': 'G'},   # A379G
        {'CHROM': 'CYP51', 'POS': 1413, 'REF': 'A', 'ALT': 'G'},   # I381V
        {'CHROM': 'CYP51', 'POS': 1648, 'REF': 'A', 'ALT': 'C'},   # Y459S
        {'CHROM': 'CYP51', 'POS': 1645, 'REF': 'ACTATGG', 'ALT': 'A'},# ΔY459/G460
        {'CHROM': 'CYP51', 'POS': 1654, 'REF': 'A', 'ALT': 'C'},   # Y461S
        {'CHROM': 'CYP51', 'POS': 1653, 'REF': 'T', 'ALT': 'C'},   # Y461H
        {'CHROM': 'CYP51', 'POS': 1811, 'REF': 'C', 'ALT': 'A'},   # N513K_AAA
        {'CHROM': 'CYP51', 'POS': 1811, 'REF': 'C', 'ALT': 'G'},   # N513K_AAG
        {'CHROM': 'CYP51', 'POS': 1843, 'REF': 'G', 'ALT': 'C'},    # S524T
        {'CHROM': 'CYP51', 'POS': 1500, 'REF': 'G', 'ALT': 'A'}    # A410T

    ],

    'CytB': [
        {'CHROM': 'CytB', 'POS': 110, 'REF': 'G', 'ALT': 'T'},     # G37V
        {'CHROM': 'CytB', 'POS': 428, 'REF': 'G', 'ALT': 'C'}      # G143A
    ],

    'SdhB': [
        {'CHROM': 'SdhB', 'POS': 788, 'REF': 'A', 'ALT': 'C'},     # N225T
        {'CHROM': 'SdhB', 'POS': 788, 'REF': 'A', 'ALT': 'T'},     # N225I
        {'CHROM': 'SdhB', 'POS': 908, 'REF': 'G', 'ALT': 'C'},     # R265P
        {'CHROM': 'SdhB', 'POS': 914, 'REF': 'A', 'ALT': 'T'},     # H267L
        {'CHROM': 'SdhB', 'POS': 917, 'REF': 'C', 'ALT': 'T'},     # T268I
        {'CHROM': 'SdhB', 'POS': 916, 'REF': 'A', 'ALT': 'G'}      # T268A
    ],

    'SdhC': [
        {'CHROM': 'SdhC', 'POS': 365, 'REF': 'C', 'ALT': 'A'},     # T79N
        {'CHROM': 'SdhC', 'POS': 368, 'REF': 'G', 'ALT': 'C'},     # W80S
        {'CHROM': 'SdhC', 'POS': 379, 'REF': 'G', 'ALT': 'T'},     # A84F
        {'CHROM': 'SdhC', 'POS': 386, 'REF': 'A', 'ALT': 'G'},     # N86S
        {'CHROM': 'SdhC', 'POS': 387, 'REF': 'C', 'ALT': 'A'},     # N86K_AAA
        {'CHROM': 'SdhC', 'POS': 387, 'REF': 'C', 'ALT': 'G'},     # N86K_AAG
        {'CHROM': 'SdhC', 'POS': 508, 'REF': 'C', 'ALT': 'G'},     # P127A
        {'CHROM': 'SdhC', 'POS': 582, 'REF': 'G', 'ALT': 'C'},     # R151S_AGC
        {'CHROM': 'SdhC', 'POS': 582, 'REF': 'G', 'ALT': 'T'},     # R151S_AGT
        {'CHROM': 'SdhC', 'POS': 581, 'REF': 'G', 'ALT': 'C'},     # R151T
        {'CHROM': 'SdhC', 'POS': 581, 'REF': 'G', 'ALT': 'T'},     # R151M
        {'CHROM': 'SdhC', 'POS': 580, 'REF': 'A', 'ALT': 'G'},     # R151G
        {'CHROM': 'SdhC', 'POS': 584, 'REF': 'A', 'ALT': 'G'},     # H152R
        {'CHROM': 'SdhC', 'POS': 611, 'REF': 'T', 'ALT': 'G'},     # I161S
        {'CHROM': 'SdhC', 'POS': 625, 'REF': 'G', 'ALT': 'A'},     # V166M
        {'CHROM': 'SdhC', 'POS': 632, 'REF': 'C', 'ALT': 'G'}      # T168R
    ],

    'SdhD': [
        {'CHROM': 'SdhD', 'POS': 139, 'REF': 'C', 'ALT': 'T'},     # R47W
        {'CHROM': 'SdhD', 'POS': 148, 'REF': 'A', 'ALT': 'T'},     # I50F
        {'CHROM': 'SdhD', 'POS': 148, 'REF': 'A', 'ALT': 'C'},     # I50L
        {'CHROM': 'SdhD', 'POS': 399, 'REF': 'A', 'ALT': 'G'},     # M114V
        {'CHROM': 'SdhD', 'POS': 497, 'REF': 'A', 'ALT': 'G'},     # D129G
        {'CHROM': 'SdhD', 'POS': 498, 'REF': 'C', 'ALT': 'A'},     # D129E_GAA
        {'CHROM': 'SdhD', 'POS': 498, 'REF': 'C', 'ALT': 'G'}      # D129E_GAG
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
    max_reads=max_reads,
)
workflow.run()