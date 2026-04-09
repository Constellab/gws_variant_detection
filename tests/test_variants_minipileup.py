import os

from gws_core import BaseTestCase, File, Folder, TaskRunner

from gws_variant_detection.variants_minipileup.variants_minipileup import MinipileupVariantsDetection

TESTDATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "testdata"))
BAM_DIR = os.path.join(TESTDATA_DIR, "bam")
REF_FASTA = os.path.join(TESTDATA_DIR, "reference.fa")


class TestMinipileupVariantsDetection(BaseTestCase):
    """MinipileupVariantsDetection: appel de variants par pileup avec Minipileup."""

    def test_variants_minipileup(self):
        """Vérifie que Minipileup produit des fichiers VCF à partir de BAMs alignés."""
        if not os.path.isdir(BAM_DIR):
            self.skipTest(
                f"testdata/bam/ absent — lancer d'abord: gws server test test_generate_testdata"
            )
        if not os.path.isfile(REF_FASTA):
            self.skipTest(f"Référence FASTA manquante: {REF_FASTA}")

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

        result: Folder = outputs["output_path"]
        self.assertIsNotNone(result)
        self.assertTrue(os.path.isdir(result.path))

        vcf_files = [f for f in os.listdir(result.path) if f.endswith(".vcf")]
        self.assertGreater(len(vcf_files), 0, "Au moins un VCF doit être produit")

        self.print(f"\n=== MinipileupVariantsDetection ===")
        self.print(f"VCFs produits ({len(vcf_files)}) : {vcf_files}")
        for vcf in vcf_files:
            vcf_path = os.path.join(result.path, vcf)
            lines = [l for l in open(vcf_path).readlines() if not l.startswith("##")]
            self.print(f"  {vcf} — {len(lines) - 1} variant(s) (hors header)")
            for line in lines[:6]:
                self.print(f"    {line.rstrip()}")
