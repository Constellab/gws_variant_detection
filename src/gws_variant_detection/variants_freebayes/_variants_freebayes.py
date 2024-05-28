import os
import subprocess
from multiprocessing import Pool
import signal
import sys

class VariantCaller:
    def __init__(self, bam_folder, ref_genome, bed_file, output_folder, ploidy, num_processes, timeout_hours):
        self.bam_folder = bam_folder
        self.ref_genome = ref_genome
        self.bed_file = bed_file
        self.output_folder = output_folder
        self.ploidy = ploidy
        self.num_processes = num_processes
        self.timeout_seconds = timeout_hours * 3600  # Convert hours to seconds

    def index_reference(self):
        subprocess.run(["minimap2", "-d", f"{self.ref_genome}.mmi", self.ref_genome], check=True)
        subprocess.run(["samtools", "faidx", self.ref_genome], check=True)

    def call_variant(self, bam_file):
        vcf_output_file = os.path.join(self.output_folder, f"{os.path.splitext(os.path.basename(bam_file))[0]}_variants_raw.vcf")
        freebayes_cmd = [
            'freebayes',
            '-f', self.ref_genome,
            '-t', self.bed_file,
            '--min-coverage', '1',
            '--max-complex-gap', '600',
            '-b', bam_file,
            '--vcf', vcf_output_file,
            '--ploidy', str(self.ploidy),
            '--min-alternate-fraction', '0'
        ]

        # Execute the command with timeout handling
        with subprocess.Popen(freebayes_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as process:
            try:
                stdout, stderr = process.communicate(timeout=self.timeout_seconds)
            except subprocess.TimeoutExpired:
                process.kill()  # Send SIGKILL to the process group
                stdout, stderr = process.communicate()
                print(f"Processing of {bam_file} reached timeout and was terminated.")

    def call_variants(self):
        self.index_reference()

        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        bam_files = [os.path.join(self.bam_folder, f) for f in os.listdir(self.bam_folder) if f.endswith('.sorted.bam')]

        with Pool(self.num_processes) as pool:
            pool.map(self.call_variant, bam_files)

bam_folder = sys.argv[1]
ref_genome = sys.argv[2]
bed_file = sys.argv[3]
output_folder = sys.argv[4]
ploidy = int(sys.argv[5])
num_processes = int(sys.argv[6])
timeout_hours = 4  # Set the timeout limit to 5 hours

variant_caller = VariantCaller(bam_folder, ref_genome, bed_file, output_folder, ploidy, num_processes, timeout_hours)
variant_caller.call_variants()
