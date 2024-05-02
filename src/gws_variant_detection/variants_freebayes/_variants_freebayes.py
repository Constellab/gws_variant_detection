import sys
import os
import subprocess

class VariantCaller:
    def __init__(self, bam_folder, ref_genome, output_folder, ploidy=2, detection_threshold=0.05):
        self.bam_folder = bam_folder
        self.ref_genome = ref_genome
        self.output_folder = output_folder
        self.ploidy = ploidy
        self.detection_threshold = detection_threshold

    def index_reference(self):
        # Index the reference genome using minimap2
        subprocess.run(["minimap2", "-d", f"{self.ref_genome}.mmi", self.ref_genome])

        # Index the reference genome using samtools faidx
        subprocess.run(["samtools", "faidx", self.ref_genome])

    def call_variants(self):
        # Index the reference genome
        self.index_reference()

        # Create the output folder if it doesn't exist
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        # Generate a list of BAM files
        bam_list_file = os.path.join(self.output_folder, "BAM.list")
        with open(bam_list_file, "w") as f:
            for bam_file in os.listdir(self.bam_folder):
                if bam_file.endswith('.dedup.bam'):
                    f.write(os.path.join(self.bam_folder, bam_file) + '\n')

        # Run FreeBayes with the detection threshold
        vcf_output_file = os.path.join(self.output_folder, "variants_raw.vcf")
        freebayes_cmd = f'freebayes -f {self.ref_genome} -L {bam_list_file} --vcf {vcf_output_file} --ploidy {self.ploidy} --min-alternate-fraction {self.detection_threshold}'
        subprocess.run(freebayes_cmd, shell=True)


# Example usage
bam_folder = sys.argv[1]
ref_genome = sys.argv[2]
output_folder = sys.argv[3]
ploidy = int(sys.argv[4])
detection_threshold = float(sys.argv[5])

variant_caller = VariantCaller(bam_folder, ref_genome, output_folder, ploidy, detection_threshold)
variant_caller.call_variants()
