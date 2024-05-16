import os
import subprocess
import sys

class MinipileupCaller:
    def __init__(self, minipileup_file_path , bam_folder, ref_genome, output_folder):
        self.minipileup_file_path = minipileup_file_path
        self.bam_folder = bam_folder
        self.ref_genome = ref_genome
        self.output_folder = output_folder

    def index_reference(self):
        # Index the reference genome using samtools faidx, if required
        subprocess.run(["samtools", "faidx", self.ref_genome])

    def call_variants(self):
        # Index the reference genome
        self.index_reference()

        # Create the output folder if it doesn't exist
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        # Iterate over each BAM file in the folder
        for bam_file in os.listdir(self.bam_folder):
            if bam_file.endswith('.sorted.bam'):
                bam_path = os.path.join(self.bam_folder, bam_file)
                output_vcf = os.path.join(self.output_folder, f"{os.path.splitext(bam_file)[0]}_variants.vcf")
                # Run Minipileup for each BAM file
                minipileup_cmd = f'{self.minipileup_file_path} -vc -f {self.ref_genome} {bam_path} > {output_vcf}'
                subprocess.run(minipileup_cmd, shell=True)

# Example usage
minipileup_file_path = sys.argv[1]
bam_folder = sys.argv[2]
ref_genome = sys.argv[3]
output_folder = sys.argv[4]

minipileup_caller = MinipileupCaller(minipileup_file_path , bam_folder, ref_genome, output_folder)
minipileup_caller.call_variants()
