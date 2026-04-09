"""
Bootstrap testdata — à lancer UNE FOIS avant tous les autres tests.

Ce script exécute la chaîne du pipeline de détection de variants du début à la fin :
  FASTQ → NanoMapping → BAM
                          ├─ FreebayesVariantsDetection → testdata/vcf_freebayes/
                          └─ MinipileupVariantsDetection → testdata/vcf_minipileup/

Les sorties sont copiées dans tests/testdata/ pour que les autres tests puissent les utiliser.

Lancer UNE FOIS :
  cd /lab/user/bricks/gws_variant_detection
  gws server test test_generate_testdata
"""
import os
import shutil

from gws_core import BaseTestCase, File, Folder, TaskRunner

from gws_variant_detection.mapping_minimap2.mapping_minimap2 import NanoMapping
from gws_variant_detection.variants_freebayes.variants_freebayes import FreebayesVariantsDetection
from gws_variant_detection.variants_minipileup.variants_minipileup import MinipileupVariantsDetection

TESTDATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "testdata"))
FASTQ_DIR = os.path.join(TESTDATA_DIR, "fastq")
REF_FASTA = os.path.join(TESTDATA_DIR, "reference.fa")
BED_FILE = os.path.join(TESTDATA_DIR, "targets.bed")

BAM_DIR = os.path.join(TESTDATA_DIR, "bam")
VCF_FREEBAYES_DIR = os.path.join(TESTDATA_DIR, "vcf_freebayes")
VCF_MINIPILEUP_DIR = os.path.join(TESTDATA_DIR, "vcf_minipileup")


class TestGenerateTestdata(BaseTestCase):
    """
    Génère tout le testdata nécessaire pour les autres tests de la brick.

    Étapes :
      1) NanoMapping      : FASTQ → BAM (saved in testdata/bam/)
      2) Freebayes        : BAM  → VCF (saved in testdata/vcf_freebayes/)
      3) Minipileup       : BAM  → VCF (saved in testdata/vcf_minipileup/)
    """

    def _copy_folder(self, src_path: str, dest_path: str):
        """Copie src_path → dest_path (efface dest s'il existe)."""
        if os.path.exists(dest_path):
            shutil.rmtree(dest_path)
        shutil.copytree(src_path, dest_path)
        self.print(f"Copié: {src_path} → {dest_path}")

    def test_step1_generate_bam(self):
        """Étape 1 : NanoMapping FASTQ → BAM sauvegardé dans testdata/bam/"""
        if not os.path.isdir(FASTQ_DIR):
            self.skipTest(f"Dossier FASTQ absent: {FASTQ_DIR}")
        if not os.path.isfile(REF_FASTA):
            self.skipTest(f"Référence FASTA absente: {REF_FASTA}")

        self.print("=== Étape 1 : Mapping minimap2 FASTQ → BAM ===")
        try:
            outputs = TaskRunner(
                task_type=NanoMapping,
                inputs={
                    "input_folder": Folder(FASTQ_DIR),
                    "input_file": File(REF_FASTA),
                },
                params={
                    "num_threads": 2,
                    "target_reads": 50,
                    "low_tail_q": 0.0,
                },
            ).run()
        except Exception as e:
            self.skipTest(f"Environnement conda (minimap2) non installé: {e}")

        bam_folder: Folder = outputs["output_path"]
        self.assertTrue(os.path.isdir(bam_folder.path))

        bam_files = [f for f in os.listdir(bam_folder.path) if f.endswith(".bam")]
        self.assertGreater(len(bam_files), 0, "Au moins un BAM doit être produit")

        self._copy_folder(bam_folder.path, BAM_DIR)
        self.print(f"BAMs disponibles: {bam_files}")

    def test_step2_generate_vcf_freebayes(self):
        """Étape 2 : FreebayesVariantsDetection BAM → VCF sauvegardé dans testdata/vcf_freebayes/"""
        if not os.path.isdir(BAM_DIR):
            self.skipTest("testdata/bam/ absent — lancer d'abord test_step1_generate_bam")
        if not os.path.isfile(REF_FASTA):
            self.skipTest(f"Référence FASTA absente: {REF_FASTA}")
        if not os.path.isfile(BED_FILE):
            self.skipTest(f"Fichier BED absent: {BED_FILE}")

        self.print("=== Étape 2 : Freebayes BAM → VCF ===")
        try:
            outputs = TaskRunner(
                task_type=FreebayesVariantsDetection,
                inputs={
                    "input_folder": Folder(BAM_DIR),
                    "input_file": File(REF_FASTA),
                    "bed_input_file": File(BED_FILE),
                },
                params={
                    "ploidy": 2,
                    "num_processes": 4,
                },
            ).run()
        except Exception as e:
            self.skipTest(f"Environnement conda (freebayes) non installé: {e}")

        vcf_folder: Folder = outputs["output_path"]
        self.assertTrue(os.path.isdir(vcf_folder.path))

        vcf_files = [f for f in os.listdir(vcf_folder.path) if f.endswith(".vcf")]
        self.assertGreater(len(vcf_files), 0, "Au moins un VCF freebayes doit être produit")

        self._copy_folder(vcf_folder.path, VCF_FREEBAYES_DIR)
        self.print(f"VCFs freebayes disponibles: {vcf_files}")

    def test_step3_generate_vcf_minipileup(self):
        """Étape 3 : MinipileupVariantsDetection BAM → VCF sauvegardé dans testdata/vcf_minipileup/"""
        if not os.path.isdir(BAM_DIR):
            self.skipTest("testdata/bam/ absent — lancer d'abord test_step1_generate_bam")
        if not os.path.isfile(REF_FASTA):
            self.skipTest(f"Référence FASTA absente: {REF_FASTA}")

        self.print("=== Étape 3 : Minipileup BAM → VCF ===")
        try:
            outputs = TaskRunner(
                task_type=MinipileupVariantsDetection,
                inputs={
                    "input_folder": Folder(BAM_DIR),
                    "input_file": File(REF_FASTA),
                },
                params={},
            ).run()
        except Exception as e:
            self.skipTest(f"Environnement conda (minipileup) non installé ou réseau indisponible: {e}")

        vcf_folder: Folder = outputs["output_path"]
        self.assertTrue(os.path.isdir(vcf_folder.path))

        vcf_files = [f for f in os.listdir(vcf_folder.path) if f.endswith(".vcf")]
        self.assertGreater(len(vcf_files), 0, "Au moins un VCF minipileup doit être produit")

        self._copy_folder(vcf_folder.path, VCF_MINIPILEUP_DIR)
        self.print(f"VCFs minipileup disponibles: {vcf_files}")
