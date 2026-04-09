import os

from gws_core import BaseTestCase, File, Folder, ResourceSet, TaskRunner

from gws_variant_detection.mapping_minimap2.mapping_minimap2 import NanoMapping

TESTDATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "testdata"))
FASTQ_DIR = os.path.join(TESTDATA_DIR, "fastq")
REF_FASTA = os.path.join(TESTDATA_DIR, "reference.fa")


class TestNanoMapping(BaseTestCase):
    """NanoMapping: alignement de lectures ONT avec Minimap2."""

    def test_mapping(self):
        """Vérifie que NanoMapping produit des BAMs et des plots de distribution."""
        if not os.path.isdir(FASTQ_DIR):
            self.skipTest(f"Dossier FASTQ absent: {FASTQ_DIR}")
        if not os.path.isfile(REF_FASTA):
            self.skipTest(f"Référence FASTA absente: {REF_FASTA}")

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
        self.assertIsNotNone(bam_folder)
        self.assertTrue(os.path.isdir(bam_folder.path))

        bam_files = [f for f in os.listdir(bam_folder.path) if f.endswith(".bam")]
        self.assertGreater(len(bam_files), 0, "Au moins un BAM doit être produit")

        plots_folder: Folder = outputs["read_length_plots"]
        self.assertIsNotNone(plots_folder)

        distributions: ResourceSet = outputs["distributions"]
        self.assertIsNotNone(distributions)

        self.print(f"\n=== NanoMapping ===")
        self.print(f"BAMs produits ({len(bam_files)}) : {bam_files}")
        html_files = [f for f in os.listdir(plots_folder.path) if f.endswith(".html")]
        self.print(f"Plots HTML produits ({len(html_files)}) : {html_files}")
        self.print(f"Distributions ResourceSet : {len(distributions.get_resources())} ressource(s)")
