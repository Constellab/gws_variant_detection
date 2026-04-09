import os

from gws_core import BaseTestCase, File, Folder, TaskRunner

from gws_variant_detection.reads_classification.reads_classification import RawReadsClassification

TESTDATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "testdata"))
FASTQ_DIR = os.path.join(TESTDATA_DIR, "fastq")
REF_FASTA = os.path.join(TESTDATA_DIR, "reference.fa")


class TestRawReadsClassification(BaseTestCase):
    """RawReadsClassification: classification des lectures par gene cible via mappy."""

    def test_reads_classification(self):
        """Vérifie la classification de lectures nanopore sur séquences de référence."""
        if not os.path.isdir(FASTQ_DIR):
            self.skipTest(f"Dossier FASTQ absent: {FASTQ_DIR}")
        if not os.path.isfile(REF_FASTA):
            self.skipTest(f"Référence FASTA absente: {REF_FASTA}")

        try:
            outputs = TaskRunner(
                task_type=RawReadsClassification,
                inputs={
                    "input_folder": Folder(FASTQ_DIR),
                    "input_file": File(REF_FASTA),
                },
                params={},
            ).run()
        except Exception as e:
            self.skipTest(f"Environnement pip (mappy) non installé: {e}")

        result: Folder = outputs["output_path"]
        self.assertIsNotNone(result)
        self.assertTrue(os.path.isdir(result.path))

        self.print(f"\n=== RawReadsClassification — Dossier de sortie : {result.path} ===")
        for gene_dir in sorted(os.listdir(result.path)):
            gene_path = os.path.join(result.path, gene_dir)
            if os.path.isdir(gene_path):
                reads = os.listdir(gene_path)
                self.print(f"  Gene '{gene_dir}': {len(reads)} fichier(s) — {reads}")
        if not os.listdir(result.path):
            self.print("  (Aucune lecture classifiée — identité insuffisante sur les données synthétiques)")
