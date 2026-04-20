import os

from gws_core import BaseTestCase, File, Folder, TaskRunner

from gws_variant_detection.variants_freebayes.variants_freebayes import FreebayesVariantsDetection

TESTDATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "testdata"))
BAM_DIR = os.path.join(TESTDATA_DIR, "bam")
REF_FASTA = os.path.join(TESTDATA_DIR, "reference.fa")
BED_FILE = os.path.join(TESTDATA_DIR, "targets.bed")


class TestFreebayesVariantsDetection(BaseTestCase):
    """FreebayesVariantsDetection: appel de variants bayésien avec Freebayes."""

    def test_variants_freebayes(self):
        """Vérifie que Freebayes produit des fichiers VCF à partir de BAMs alignés."""
        if not os.path.isdir(BAM_DIR):
            self.skipTest(
                f"testdata/bam/ absent — lancer d'abord: gws server test test_generate_testdata"
            )
        if not os.path.isfile(REF_FASTA):
            self.skipTest(f"Référence FASTA manquante: {REF_FASTA}")
        if not os.path.isfile(BED_FILE):
            self.skipTest(f"Fichier BED manquant: {BED_FILE}")

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

        result: Folder = outputs["output_path"]
        self.assertIsNotNone(result)
        self.assertTrue(os.path.isdir(result.path))

        vcf_files = [f for f in os.listdir(result.path) if f.endswith(".vcf")]
        self.assertGreater(len(vcf_files), 0, "Au moins un VCF doit être produit")

        self.print(f"\n=== FreebayesVariantsDetection ===")
        self.print(f"VCFs produits ({len(vcf_files)}) : {vcf_files}")
        for vcf in vcf_files:
            vcf_path = os.path.join(result.path, vcf)
            lines = [l for l in open(vcf_path).readlines() if not l.startswith("##")]
            self.print(f"  {vcf} — {len(lines) - 1} variant(s) (hors header)")
            for line in lines[:6]:
                self.print(f"    {line.rstrip()}")
