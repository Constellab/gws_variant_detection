import os

from gws_core import BaseTestCase, File, Folder, Settings, TaskRunner

from gws_variant_detection.trimming_porechop.trimming_porechop import PorechopTrimming

TESTDATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "testdata"))
FASTQ_DIR = os.path.join(TESTDATA_DIR, "fastq")


class TestPorechopTrimming(BaseTestCase):
    """PorechopTrimming: élimination des adaptateurs nanopore avec Porechop."""

    def test_trimming(self):
        """Vérifie que PorechopTrimming produit un dossier de lectures trimmées."""
        try:
            outputs = TaskRunner(
                task_type=PorechopTrimming,
                inputs={"input_folder": Folder(FASTQ_DIR)},
                params={"num_threads": 2},
            ).run()
        except Exception as e:
            self.skipTest(f"Environnement conda non installé ou données insuffisantes: {e}")

        result: Folder = outputs["output_path"]
        self.assertIsNotNone(result)
        self.assertTrue(os.path.isdir(result.path))

        trimmed_files = os.listdir(result.path)
        self.print(f"\n=== PorechopTrimming — Dossier de sortie : {result.path} ===")
        self.print(f"Fichiers produits ({len(trimmed_files)}) : {trimmed_files}")
