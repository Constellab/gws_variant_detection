import os
import subprocess
from multiprocessing import Pool
import sys

class VariantCaller:
    def __init__(self, bam_folder, ref_genome, output_folder, ploidy, num_processes, bed_file):
        self.bam_folder = bam_folder
        self.ref_genome = ref_genome
        self.output_folder = output_folder
        self.ploidy = ploidy
        self.num_processes = num_processes
        self.bed_file = bed_file

    def index_reference(self):
        # Index the reference genome using minimap2
        subprocess.run(["minimap2", "-d", f"{self.ref_genome}.mmi", self.ref_genome])

        # Index the reference genome using samtools faidx
        subprocess.run(["samtools", "faidx", self.ref_genome])

    def call_variant(self, bam_file):
        vcf_output_file = os.path.join(self.output_folder, f"{os.path.splitext(os.path.basename(bam_file))[0]}_variants_raw.vcf")
        freebayes_cmd = f'freebayes -f {self.ref_genome} -t "{self.bed_file}" --min-coverage 1 --max-complex-gap 600 -b {bam_file} --vcf {vcf_output_file} --ploidy {self.ploidy} --min-alternate-fraction 0'
        subprocess.run(freebayes_cmd, shell=True)

    def call_variants(self):
        # Index the reference genome
        self.index_reference()

        # Create the output folder if it doesn't exist
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        # Generate a list of BAM files
        bam_files = [os.path.join(self.bam_folder, f) for f in os.listdir(self.bam_folder) if f.endswith('.sorted.bam')]

        # Call variants in parallel
        with Pool(self.num_processes) as pool:
            pool.map(self.call_variant, bam_files)


# Example usage
bam_folder = sys.argv[1]
ref_genome = sys.argv[2]
bed_file = sys.argv[3]
output_folder = sys.argv[4]
ploidy = int(sys.argv[5])
num_processes = int(sys.argv[6])

variant_caller = VariantCaller(bam_folder, ref_genome, output_folder, ploidy, num_processes, bed_file)
variant_caller.call_variants()
