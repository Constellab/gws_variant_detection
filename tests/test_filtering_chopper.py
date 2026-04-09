import os

from gws_core import BaseTestCase, File, Folder, ResourceSet, TaskRunner
from gws_core import Table

from gws_variant_detection.filtering_chopper.filtering_chopper import ReadsQualityFiltering

TESTDATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "testdata"))
FASTQ_DIR = os.path.join(TESTDATA_DIR, "fastq")
REF_FASTA = os.path.join(TESTDATA_DIR, "reference.fa")


class TestReadsQualityFiltering(BaseTestCase):
    """ReadsQualityFiltering: filtrage Phred + identité/couverture avec chopper et mappy."""

    def test_filtering_chopper(self):
        """Vérifie que le filtrage produit des lectures filtrées et des statistiques."""
        if not os.path.isdir(FASTQ_DIR):
            self.skipTest(f"Dossier FASTQ absent: {FASTQ_DIR}")
        if not os.path.isfile(REF_FASTA):
            self.skipTest(f"Référence FASTA absente: {REF_FASTA}")

        try:
            outputs = TaskRunner(
                task_type=ReadsQualityFiltering,
                inputs={
                    "input_folder": Folder(FASTQ_DIR),
                    "reference_fasta": File(REF_FASTA),
                },
                params={
                    "num_threads": 2,
                    "qualified_quality_phred": 5,
                    "min_identity": 50,
                    "min_coverage": 50,
                },
            ).run()
        except Exception as e:
            self.skipTest(f"Environnement conda (chopper/mappy) non installé: {e}")

        filtered_folder: Folder = outputs["final_filtered_reads"]
        self.assertIsNotNone(filtered_folder)
        self.assertTrue(os.path.isdir(filtered_folder.path))

        summary_table: Table = outputs["number_reads_filtered"]
        self.assertIsNotNone(summary_table)
        df = summary_table.to_dataframe()
        self.assertGreater(df.shape[1], 0)

        plots_raw: Folder = outputs["plots_raw"]
        self.assertIsNotNone(plots_raw)

        rs_raw: ResourceSet = outputs["plots_raw_rs"]
        self.assertIsNotNone(rs_raw)

        filtered_files = os.listdir(filtered_folder.path)
        self.print(f"\n=== ReadsQualityFiltering ===")
        self.print(f"Lectures filtrées ({len(filtered_files)}) : {filtered_files}")
        self.print(f"\nTable résumé du filtrage ({df.shape[0]} lignes x {df.shape[1]} colonnes) :")
        self.print(df.to_string())
        raw_html = [f for f in os.listdir(plots_raw.path) if f.endswith(".html")]
        self.print(f"\nPlots RAW produits ({len(raw_html)}) : {raw_html}")
        self.print(f"ResourceSet RAW : {len(rs_raw.get_resources())} ressource(s)")
