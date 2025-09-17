#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import json
import re
import ast
import shutil
import subprocess
from multiprocessing import Pool

import pandas as pd
import pysam


class VCFWorkflowFiltering:
    """
    Septoria variant workflow driven by the same JSON schema.

    See _variant_norm_filtering.py docstring for the JSON structure.
    """

    def __init__(
        self,
        bam_folder: str,
        ref_genome: str,
        gene_name: str,
        config_json_path: str,
        minipileup_path: str,
        output_folder: str,
        threads: int,
        max_reads: int,
    ):
        self.bam_folder = bam_folder
        self.ref_genome = ref_genome
        self.gene_name = gene_name
        self.minipileup_path = minipileup_path
        self.output_folder = output_folder
        self.threads = threads
        self.max_reads = max_reads

        # Folders
        self._mk(self.output_folder)
        self.per_read_bam_folder = os.path.join(self.output_folder, "per_read_bams")
        self.per_read_vcf_folder = os.path.join(self.output_folder, "per_read_vcfs")
        self.normalized_vcf_folder = os.path.join(self.output_folder, "normalized_vcfs")
        self._mk(self.per_read_bam_folder)
        self._mk(self.per_read_vcf_folder)
        self._mk(self.normalized_vcf_folder)

        # Load config JSON and select gene block (CYP51, CytB, SdhB, SdhC, SdhD, etc.)
        with open(config_json_path, "r") as f:
            cfg_all = json.load(f)
        if "genes" not in cfg_all or self.gene_name not in cfg_all["genes"]:
            raise ValueError(f"Gene '{self.gene_name}' not found in {config_json_path}")
        self.gene_cfg = cfg_all["genes"][self.gene_name]

        # Map (chrom,pos,ref,alt) -> [name,...] for simple sites
        self.site_to_names = {}
        for section in ("substitution_simple", "insertion_simple", "deletion_simple"):
            for entry in self.gene_cfg.get(section, []):
                key = (self.gene_name, int(entry["pos"]), entry["ref"], entry["alt"])
                self.site_to_names.setdefault(key, []).append(entry["name"])

        # Complex patterns
        self.insertion_signatures = self.gene_cfg.get("insertion_signature", [])
        self.substitution_composites = self.gene_cfg.get("substitution_composite", [])

        # Known names
        self.known_variant_names = set()
        for section in ("substitution_simple", "insertion_simple", "deletion_simple"):
            for e in self.gene_cfg.get(section, []):
                self.known_variant_names.add(e["name"])
        for e in self.insertion_signatures:
            self.known_variant_names.add(e["name"])
        for e in self.substitution_composites:
            self.known_variant_names.add(e["name"])

        # Book-keeping
        self.barcode_list = self._get_barcodes()
        self.read_variant_records = []
        self.total_reads_per_barcode = {}

        # Output files
        self.read_variants_path = os.path.join(self.output_folder, "read_variants_table.csv")
        if not os.path.exists(self.read_variants_path):
            pd.DataFrame(columns=["Read", "Variants"]).to_csv(self.read_variants_path, index=False)
        self.filtered_reads_count_path = os.path.join(self.output_folder, "filtered_reads_count.csv")

    # -------------- utils --------------

    def _mk(self, p: str):
        os.makedirs(p, exist_ok=True)

    def _get_barcodes(self):
        barcodes = set()
        for bam_file in os.listdir(self.bam_folder):
            if bam_file.endswith(".bam"):
                barcode = os.path.splitext(bam_file)[0].split("-")[0]
                barcodes.add(barcode)
        return sorted(list(barcodes))

    # -------------- stage 1 --------------

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

    def _save_filtered_reads_count(self, barcode, count):
        df = pd.DataFrame([[barcode, count]], columns=["Barcode", "Filtered_Read_Count"])
        if os.path.exists(self.filtered_reads_count_path):
            df.to_csv(self.filtered_reads_count_path, mode="a", header=False, index=False)
        else:
            df.to_csv(self.filtered_reads_count_path, index=False)

    def split_bam_by_read(self):
        bam_files = [f for f in os.listdir(self.bam_folder) if f.endswith(".bam")]
        bam_paths = [os.path.join(self.bam_folder, f) for f in bam_files]
        with Pool(processes=self.threads) as pool:
            pool.map(self._process_bam_file, bam_paths)

    def _process_bam_file(self, bam_path):
        bam_file = os.path.basename(bam_path)
        bam_basename = os.path.splitext(bam_file)[0]
        barcode = bam_basename.split("-")[0]
        barcode_folder = os.path.join(self.per_read_bam_folder, barcode)
        self._mk(barcode_folder)
        try:
            filtered_reads, unique_count = self.filter_unique_reads(bam_path)
            template = pysam.AlignmentFile(bam_path, "rb")
            for idx, read in enumerate(filtered_reads, start=1):
                per_read_bam_path = os.path.join(barcode_folder, f"{bam_basename}_read{idx}.bam")
                with pysam.AlignmentFile(per_read_bam_path, "wb", template=template) as out_bam:
                    out_bam.write(read)
                pysam.index(per_read_bam_path)
            template.close()
            self._save_filtered_reads_count(barcode, unique_count)
        except Exception as e:
            print(f"Erreur lors du traitement de {bam_file} : {e}")

    # -------------- stage 2 --------------

    def call_variants_per_read(self):
        bam_folders = [os.path.join(self.per_read_bam_folder, bc) for bc in self.barcode_list]
        with Pool(processes=self.threads) as pool:
            pool.map(self._process_bam_folder, bam_folders)

    def _process_bam_folder(self, barcode_bam_folder):
        barcode = os.path.basename(barcode_bam_folder)
        barcode_vcf_folder = os.path.join(self.per_read_vcf_folder, barcode)
        self._mk(barcode_vcf_folder)
        bam_files = [f for f in os.listdir(barcode_bam_folder) if f.endswith(".bam")]
        for bam_file in bam_files:
            bam_path = os.path.join(barcode_bam_folder, bam_file)
            vcf_file = bam_file.replace(".bam", ".vcf")
            vcf_path = os.path.join(barcode_vcf_folder, vcf_file)
            cmd = f"{self.minipileup_path} -vc -f {self.ref_genome} {bam_path}"
            try:
                with open(vcf_path, "w") as vcf_output:
                    subprocess.run(cmd, shell=True, check=True, stdout=vcf_output, stderr=subprocess.PIPE)
            except subprocess.CalledProcessError as e:
                print(f"Minipileup a échoué pour {bam_file} : {e.stderr.decode().strip()}")

    # -------------- stage 3 --------------

    def normalize_vcfs(self):
        vcf_folders = [os.path.join(self.per_read_vcf_folder, bc) for bc in self.barcode_list]
        with Pool(processes=self.threads) as pool:
            pool.map(self._process_vcf_folder, vcf_folders)

    def _process_vcf_folder(self, barcode_vcf_folder):
        barcode = os.path.basename(barcode_vcf_folder)
        norm_folder = os.path.join(self.normalized_vcf_folder, barcode)
        self._mk(norm_folder)
        vcf_files = [f for f in os.listdir(barcode_vcf_folder) if f.endswith(".vcf")]
        for vcf_file in vcf_files:
            input_vcf = os.path.join(barcode_vcf_folder, vcf_file)
            step1_vcf = os.path.join(norm_folder, f"{os.path.splitext(vcf_file)[0]}_norm_step1.vcf")
            step2_vcf = os.path.join(norm_folder, f"{os.path.splitext(vcf_file)[0]}_norm_step2.vcf")
            try:
                subprocess.run(["bcftools", "norm", "-m-", "-o", step1_vcf, input_vcf], check=True)
                subprocess.run(["bcftools", "norm", "-f", self.ref_genome, "-cw", "-o", step2_vcf, step1_vcf], check=True)
            except subprocess.CalledProcessError as e:
                print(f"bcftools norm a échoué pour {vcf_file} (barcode {barcode}) : {e.stderr.decode().strip()}")

    # -------------- stage 4 --------------

    def extract_variants_from_vcf(self, vcf_path):
        present_names = set()
        present_sites = set()   # ALL observed (pos,ref,alt)
        ad_alt_by_site = {}

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

                    present_sites.add((pos, ref, alt))

                    key = (chrom, pos, ref, alt)
                    if key in self.site_to_names:
                        for nm in self.site_to_names[key]:
                            present_names.add(nm)

                    if len(cols) >= 10:
                        fmt = cols[8].split(":")
                        sample_vals = cols[9].split(":")
                        if "AD" in fmt:
                            ad_index = fmt.index("AD")
                            if ad_index < len(sample_vals):
                                parts = sample_vals[ad_index].split(",")
                                if len(parts) >= 2 and parts[1].isdigit():
                                    ad_alt_by_site[(pos, ref, alt)] = ad_alt_by_site.get((pos, ref, alt), 0) + int(parts[1])
        except Exception as e:
            print(f"Erreur lecture {vcf_path}: {e}")

        # insertion signatures
        for sig in self.insertion_signatures:
            ssum = sum(ad_alt_by_site.get((int(s["pos"]), s["ref"], s["alt"]), 0) for s in sig["sites"])
            if ssum > 0:
                present_names.add(sig["name"])

        # substitution composites
        for comp in self.substitution_composites:
            comps = [(int(s["pos"]), s["ref"], s["alt"]) for s in comp["components"]]
            if all(c in present_sites for c in comps):
                present_names.add(comp["name"])
                for c in comps:
                    k = (self.gene_name, c[0], c[1], c[2])
                    for nm in self.site_to_names.get(k, []):
                        if nm != comp["name"]:
                            present_names.discard(nm)

        return list(present_names)

    def collect_variants_and_generate_files(self):
        for barcode in self.barcode_list:
            barcode_norm_vcf_folder = os.path.join(self.normalized_vcf_folder, barcode)
            if not os.path.exists(barcode_norm_vcf_folder):
                continue
            vcf_files = [f for f in os.listdir(barcode_norm_vcf_folder) if f.endswith("_norm_step2.vcf")]
            total_processed = 0
            for vcf_file in vcf_files:
                vcf_path = os.path.join(barcode_norm_vcf_folder, vcf_file)
                read_variants = self.extract_variants_from_vcf(vcf_path)
                total_processed += 1
                read_id = os.path.splitext(vcf_file)[0].replace("_norm_step2", "")
                self.read_variant_records.append({"Read": read_id, "Variants": read_variants})
            self.total_reads_per_barcode[barcode] = total_processed

        self._save_read_variants_table()
        self._generate_known_variants_filtered()
        self._generate_unknown_variants()

    def _save_read_variants_table(self):
        df = pd.DataFrame(self.read_variant_records)
        if not df.empty:
            df["Variants"] = df["Variants"].apply(lambda vs: str(vs) if vs else "[]")
            df.to_csv(self.read_variants_path, index=False)
        else:
            print("Aucune donnée de variantes à sauvegarder dans read_variants_table.csv")

    # -------------- known / unknown / wildtype --------------

    def _calc_total_reads_per_barcode(self):
        totals = {}
        if os.path.exists(self.read_variants_path):
            df = pd.read_csv(self.read_variants_path)
            if not df.empty and "Read" in df.columns:
                df["Barcode"] = df["Read"].apply(lambda x: x.split("_")[0].split("-")[0])
                totals = df.groupby("Barcode").size().to_dict()
        return totals

    def _generate_known_variants_filtered(self):
        if not os.path.exists(self.read_variants_path):
            print("read_variants_table.csv introuvable.")
            return
        df = pd.read_csv(self.read_variants_path)
        if df.empty or "Read" not in df.columns or "Variants" not in df.columns:
            print("read_variants_table.csv vide ou colonnes manquantes.")
            return

        df["Barcode"] = df["Read"].apply(lambda x: x.split("_")[0].split("-")[0])
        df["Variants"] = df["Variants"].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else [])
        df["Known_Variants"] = df["Variants"].apply(lambda vs: [v for v in vs if v in self.known_variant_names])

        known_exploded = df.explode("Known_Variants").dropna(subset=["Known_Variants"])

        out_path = os.path.join(self.output_folder, "known_variants_filtered.tsv")
        if known_exploded.empty:
            pd.DataFrame(columns=["variant_name"]).to_csv(out_path, sep="\t", index=False)
            print("Aucun variant connu détecté – fichier vide créé.")
            return

        known_counts = (
            known_exploded.groupby(["Barcode", "Known_Variants"]).size().reset_index(name="Count")
        )

        pivot_counts = (
            known_counts.pivot_table(index="Known_Variants", columns="Barcode", values="Count", fill_value=0)
            .reset_index()
            .rename(columns={"Known_Variants": "variant_name"})
        )

        self.total_reads_per_barcode = self._calc_total_reads_per_barcode()
        for bc in self.barcode_list:
            count_col = f"{bc}_VARIANT_COUNT"
            freq_col = f"{bc}_VARIANT_FREQ"
            total = self.total_reads_per_barcode.get(bc, 1)
            if bc in pivot_counts.columns:
                pivot_counts[freq_col] = pivot_counts[bc].apply(lambda c: f"{(c / total * 100):.2f}%" if total > 0 else "0.00%")
                pivot_counts.rename(columns={bc: count_col}, inplace=True)
            else:
                pivot_counts[count_col] = 0
                pivot_counts[freq_col] = "0.00%"

        ordered_cols = ["variant_name"]
        for bc in self.barcode_list:
            ordered_cols += [f"{bc}_VARIANT_COUNT", f"{bc}_VARIANT_FREQ"]
        pivot_counts[ordered_cols].to_csv(out_path, sep="\t", index=False)

    def _get_variant_name(self, chrom, pos, ref, alt):
        key = (chrom, int(pos), ref, alt)
        names = self.site_to_names.get(key)
        return names[0] if names else None

    def _parse_minipileup_output_file(self, vcf_path):
        variants = []
        try:
            with open(vcf_path, "r") as vcf:
                for line in vcf:
                    if line.startswith("#"):
                        continue
                    cols = line.strip().split("\t")
                    if len(cols) < 10:
                        continue
                    chrom = cols[0]
                    pos = int(cols[1])
                    ref = cols[3]
                    alt = cols[4]
                    fmt = cols[8].split(":")
                    sample_vals = cols[9].split(":")

                    ad_value = "0,0"
                    if "AD" in fmt:
                        ad_index = fmt.index("AD")
                        if ad_index < len(sample_vals):
                            parts = sample_vals[ad_index].split(",")
                            if len(parts) >= 2 and all(x.isdigit() for x in parts[:2]):
                                ad_value = f"{int(parts[0])},{int(parts[1])}"

                    variant_id = f"{chrom}:{pos}:{ref}>{alt}"
                    variants.append({"variant": variant_id, "AD": ad_value})
        except Exception as e:
            print(f"Erreur lecture {vcf_path}: {e}")
        return variants

    def _generate_unknown_variants(self):
        print("Génération de unknown_variants.tsv via minipileup sur BAM complets…")
        unknown_variants_dict = {}

        bam_files = [f for f in os.listdir(self.bam_folder) if f.endswith(".bam")]
        for bam_file in bam_files:
            bam_path = os.path.join(self.bam_folder, bam_file)
            try:
                temp_vcf = os.path.join(self.output_folder, f"{bam_file}_temp.vcf")
                cmd = f"{self.minipileup_path} -vc -f {self.ref_genome} {bam_path}"
                with open(temp_vcf, "w") as out:
                    subprocess.run(cmd, shell=True, check=True, stdout=out, stderr=subprocess.PIPE)

                norm_temp_vcf = os.path.join(self.output_folder, f"{bam_file}_temp_norm.vcf")
                try:
                    subprocess.run(["bcftools", "norm", "-m-", "-o", norm_temp_vcf, temp_vcf], check=True)
                    subprocess.run(["bcftools", "norm", "-f", self.ref_genome, "-cw", "-o", temp_vcf, norm_temp_vcf], check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Erreur bcftools norm sur {temp_vcf}: {e.stderr.decode().strip()}")
                    continue

                variants = self._parse_minipileup_output_file(temp_vcf)
                barcode = bam_file.split("-")[0]

                for variant_info in variants:
                    variant = variant_info["variant"]
                    ad = variant_info["AD"]
                    try:
                        chrom, pos_str, ref_alt = variant.split(":")
                        pos = int(pos_str)
                        ref, alt = ref_alt.split(">")
                    except Exception:
                        continue

                    # known if in simple JSON site map (deletions included)
                    variant_name = self._get_variant_name(chrom, pos, ref, alt)
                    if variant_name is not None:
                        continue

                    key = (chrom, pos, ref, alt)
                    if key not in unknown_variants_dict:
                        unknown_variants_dict[key] = {"CHROM": chrom, "POS": pos, "REF": ref, "ALT": alt}

                    try:
                        rc, ac = map(int, ad.split(","))
                    except:
                        rc, ac = 0, 0
                    tot = rc + ac
                    freq = (ac / tot * 100) if tot > 0 else 0.0
                    unknown_variants_dict[key][barcode] = {"AD": ad, "Frequency": f"{freq:.2f}%", "Count": ac}

                os.remove(temp_vcf)
                os.remove(norm_temp_vcf)
            except subprocess.CalledProcessError as e:
                print(f"Minipileup échoué sur {bam_file}: {e.stderr.decode().strip()}")

        if not unknown_variants_dict:
            print("Aucune variante inconnue détectée.")
            return

        rows = []
        for key, data in unknown_variants_dict.items():
            row = [data["CHROM"], data["POS"], data["REF"], data["ALT"]]
            for bc in self.barcode_list:
                if bc in data:
                    row += [data[bc]["AD"], data[bc]["Frequency"], data[bc]["Count"]]
                else:
                    row += ["0,0", "0.00%", "0"]
            rows.append(row)

        header = ["CHROM", "POS", "REF", "ALT"]
        for bc in self.barcode_list:
            header += [f"{bc}_AD", f"{bc}_VARIANT_FREQ", f"{bc}_VARIANT_COUNT"]

        df_unknown = pd.DataFrame(rows, columns=header)
        out_unknown_path = os.path.join(self.output_folder, "unknown_variants.tsv")
        df_unknown.to_csv(out_unknown_path, sep="\t", index=False)

    def calculate_wildtype_frequencies(self):
        out_file = os.path.join(self.output_folder, "wildtype_frequencies.csv")
        results = []
        for bc in self.barcode_list:
            total_reads = self.total_reads_per_barcode.get(bc, 0)
            wt = 0
            for rec in self.read_variant_records:
                if rec["Read"].startswith(bc) and not rec["Variants"]:
                    wt += 1
            freq = (wt / total_reads * 100) if total_reads > 0 else 0.0
            results.append([bc, wt, f"{freq:.2f}%"])
        pd.DataFrame(results, columns=["Sample", "Count", "Frequency"]).to_csv(out_file, index=False)

    def clean_up(self):
        for d in [self.per_read_bam_folder, self.per_read_vcf_folder, self.normalized_vcf_folder]:
            if os.path.exists(d):
                shutil.rmtree(d)

    def run(self):
        self.split_bam_by_read()
        self.call_variants_per_read()
        self.normalize_vcfs()
        self.collect_variants_and_generate_files()
        self.calculate_wildtype_frequencies()
        self.clean_up()


if __name__ == "__main__":
    bam_folder = sys.argv[1]
    ref_genome = sys.argv[2]
    gene_name = str(sys.argv[3])
    minipileup_path = sys.argv[4]
    output_folder = sys.argv[5]
    threads = int(sys.argv[6])
    max_reads = int(sys.argv[7])
    config_json_path = sys.argv[8]

    workflow = VCFWorkflowFiltering(
        bam_folder=bam_folder,
        ref_genome=ref_genome,
        gene_name=gene_name,
        config_json_path=config_json_path,
        minipileup_path=minipileup_path,
        output_folder=output_folder,
        threads=threads,
        max_reads=max_reads,
    )
    workflow.run()
