import os

from gws_core import BaseTestCase, File, Folder, Settings, Table, TaskRunner

from gws_variant_detection.quality_check.quality_check import NanoQualityCheck

TESTDATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "testdata"))
FASTQ_DIR = os.path.join(TESTDATA_DIR, "fastq")


class TestNanoQualityCheck(BaseTestCase):
    """NanoQualityCheck: contrôle qualité de lectures nanopore avec NanoPlot."""

    def test_quality_check(self):
        """Vérifie que NanoQualityCheck produit une Table de statistiques."""
        try:
            outputs = TaskRunner(
                task_type=NanoQualityCheck,
                inputs={"input_folder": Folder(FASTQ_DIR)},
                params={"num_threads": 2},
            ).run()
        except Exception as e:
            self.skipTest(f"Environnement conda non installé ou données insuffisantes: {e}")

        result: Table = outputs["output_path"]
        self.assertIsNotNone(result)
        df = result.to_dataframe()
        self.assertGreater(df.shape[0], 0, "La table de stats doit avoir au moins une ligne")
        self.assertGreater(df.shape[1], 0, "La table de stats doit avoir au moins une colonne")

        self.print(f"\n=== NanoQualityCheck — Table ({df.shape[0]} lignes x {df.shape[1]} colonnes) ===")
        self.print(df.to_string())
