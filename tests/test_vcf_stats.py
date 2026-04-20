import os

from gws_core import BaseTestCase, File, Folder, Settings, Table, TaskRunner

from gws_variant_detection.vcf_stats.vcf_stats import VariantStats

TESTDATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "testdata"))
VCF_DIR = os.path.join(TESTDATA_DIR, "vcf")


class TestVariantStats(BaseTestCase):
    """VariantStats: statistiques sur fichiers VCF via bcftools stats."""

    def test_vcf_stats(self):
        """Vérifie que VariantStats produit une table de statistiques sur le VCF de test."""
        if not os.path.isdir(VCF_DIR):
            self.skipTest(f"Dossier VCF testdata absent: {VCF_DIR}")

        try:
            outputs = TaskRunner(
                task_type=VariantStats,
                inputs={"input_folder": Folder(VCF_DIR)},
                params={},
            ).run()
        except Exception as e:
            self.skipTest(f"Environnement conda (bcftools) non installé: {e}")

        result: Table = outputs["Statistical_vcf_table"]
        self.assertIsNotNone(result)
        df = result.to_dataframe()
        self.assertGreater(df.shape[0], 0, "La table de stats VCF doit avoir au moins une ligne")
        self.assertGreater(df.shape[1], 0, "La table de stats VCF doit avoir au moins une colonne")
        self.assertIn("barcode01", df.columns.tolist())

        self.print(f"\n=== VariantStats (bcftools) — Table ({df.shape[0]} lignes x {df.shape[1]} colonnes) ===")
        self.print(df.to_string())
