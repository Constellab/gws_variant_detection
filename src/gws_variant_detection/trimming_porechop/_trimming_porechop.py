import os
import subprocess
import sys

class PorechopTrimmer:
    def __init__(self, input_folder, output_folder, num_threads=1):
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.num_threads = num_threads

    def trim_reads(self):
        # Create the output folder if it doesn't exist
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        # Get a list of all FASTQ files in the input folder
        fastq_files = [f for f in os.listdir(self.input_folder) if f.endswith('.fastq') or f.endswith('.fastq.gz')]

        # Run Porechop for each FASTQ file
        for fastq_file in fastq_files:
            input_path = os.path.join(self.input_folder, fastq_file)
            output_path = os.path.join(self.output_folder, fastq_file)
            porechop_cmd = f'porechop -i {input_path} -o {output_path} -t {self.num_threads} --discard_middle'
            subprocess.run(porechop_cmd, shell=True)

# Example usage
input_folder = sys.argv[1]
num_threads = int(sys.argv[2])
output_folder = sys.argv[3]

porechop_trimmer = PorechopTrimmer(input_folder, output_folder, num_threads)
porechop_trimmer.trim_reads()
