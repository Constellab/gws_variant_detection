import sys
import os
import subprocess
import gzip
import mappy as mp
import pandas as pd

class FastpFiltering:
    def __init__(self, input_folder, fastp_output_folder, identity_coverage_output_folder, final_output_folder, summary_output_folder,
                 num_threads, qualified_quality_phred, reads_to_process, min_identity, min_coverage,
                 reference_genome, skip_read_threshold=100):
        self.input_folder = input_folder
        self.fastp_output_folder = fastp_output_folder
        self.identity_coverage_output_folder = identity_coverage_output_folder
        self.final_output_folder = final_output_folder
        self.summary_output_folder = summary_output_folder
        self.num_threads = num_threads
        self.qualified_quality_phred = qualified_quality_phred
        self.reads_to_process = reads_to_process
        self.min_identity = min_identity
        self.min_coverage = min_coverage
        self.reference_genome = reference_genome
        self.skip_read_threshold = skip_read_threshold
        self.aligner = mp.Aligner(reference_genome, preset="map-ont")

        # Ensure output folders exist
        os.makedirs(self.fastp_output_folder, exist_ok=True)
        os.makedirs(self.identity_coverage_output_folder, exist_ok=True)
        os.makedirs(self.final_output_folder, exist_ok=True)
        os.makedirs(self.summary_output_folder, exist_ok=True)

    def count_reads_in_fastq(self, fastq_file):
        """Counts the number of reads in a FASTQ file."""
        try:
            cmd = f'zcat {fastq_file} | wc -l'
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
            num_lines = int(result.stdout.strip())
            return num_lines // 4  # Each read takes 4 lines in a FASTQ file
        except subprocess.CalledProcessError as e:
            print(f"Error counting reads in {fastq_file}: {e}", file=sys.stderr)
            return 0

    def check_alignment_coverage(self, read_sequence):
        """Checks if a read meets the identity and coverage threshold when aligned to the reference."""
        for hit in self.aligner.map(read_sequence):
            query_length = hit.q_en - hit.q_st  # Get read length

            if query_length == 0 or hit.blen == 0:  # Avoid division by zero
                continue

            identity = (hit.mlen / hit.blen) * 100  # % of matching bases in alignment
            coverage = (hit.blen / query_length) * 100  # % of the read aligned

            if identity >= self.min_identity and coverage >= self.min_coverage:
                return identity, coverage, query_length  # Return values if valid

        return None, None, None  # Read rejected

    def filter_reads(self):
        """Filters reads based on Phred quality, identity, and coverage thresholds, and selects the longest reads."""
        summary_data = []

        fastq_files = [f for f in os.listdir(self.input_folder) if f.endswith('.fastq.gz')]

        for fastq_file in fastq_files:
            input_path = os.path.join(self.input_folder, fastq_file)
            fastp_output_path = os.path.join(self.fastp_output_folder, fastq_file)
            identity_coverage_output_path = os.path.join(self.identity_coverage_output_folder, fastq_file)
            final_output_path = os.path.join(self.final_output_folder, fastq_file)

            # Step 1: Run Fastp (Quality Filtering)
            fastp_cmd = f'fastp -i {input_path} -o {fastp_output_path} -A -G --qualified_quality_phred {self.qualified_quality_phred} --thread {self.num_threads}'
            result = subprocess.run(fastp_cmd, shell=True)
            if result.returncode != 0:
                print(f"Error with fastp on {fastq_file}", file=sys.stderr)
                continue

            # Step 2: Filter reads based on identity and coverage
            valid_reads = []
            final_reads = []

            with gzip.open(fastp_output_path, 'rt') as infile:
                fastq_lines = infile.readlines()

            with gzip.open(identity_coverage_output_path, 'wt') as outfile:
                for i in range(0, len(fastq_lines), 4):
                    read_name = fastq_lines[i].strip()
                    read_sequence = fastq_lines[i + 1].strip()
                    read_quality = fastq_lines[i + 3].strip()

                    valid_reads.append((read_name, read_sequence, read_quality))  # Store all reads after Phred filtering

                    identity, coverage, read_length = self.check_alignment_coverage(read_sequence)
                    if identity is not None and coverage is not None:
                        outfile.write(f"{read_name}\n{read_sequence}\n+\n{read_quality}\n")
                        final_reads.append((read_name, read_sequence, read_quality, identity, coverage, read_length))

            # Step 3: Sort reads by length (longest first)
            final_reads.sort(key=lambda x: x[5], reverse=True)  # Sort by read_length

            # Step 4: Create final filtered FASTQ (top 5000 longest reads)
            if len(final_reads) >= (self.reads_to_process - self.skip_read_threshold):
                with gzip.open(final_output_path, 'wt') as final_file:
                    for read in final_reads[:self.reads_to_process]:  # Take only the top longest reads_to_process
                        final_file.write(f"{read[0]}\n{read[1]}\n+\n{read[2]}\n")

            # Step 5: Create summary
            summary_data.append({
                "File": fastq_file,
                "Total_Reads_After_Phred": len(valid_reads),
                "Total_Reads_After_Identity_Coverage": len(final_reads),
                "Final_Filtered_Reads": min(len(final_reads), self.reads_to_process),
                "Min_Identity": self.min_identity,
                "Min_Coverage": self.min_coverage,
                "Phred_Threshold": self.qualified_quality_phred
            })

        # Step 6: Save summary as CSV
        summary_df = pd.DataFrame(summary_data)
        summary_csv_path = os.path.join(self.summary_output_folder, "filtering_summary.csv")
        summary_df.to_csv(summary_csv_path, index=False)
        print(f"Summary saved to {summary_csv_path}")


# Command-line arguments
input_folder = sys.argv[1]
fastp_output_folder = sys.argv[2]
identity_coverage_output_folder = sys.argv[3]
final_output_folder = sys.argv[4]
summary_output_folder = sys.argv[5]

num_threads = int(sys.argv[6])
qualified_quality_phred = int(sys.argv[7])
reads_to_process = int(sys.argv[8])
min_identity = int(sys.argv[9])
min_coverage = int(sys.argv[10])
reference_genome = sys.argv[11]
skip_read_threshold = int(sys.argv[12])

fastp_filtering = FastpFiltering(
    input_folder, fastp_output_folder, identity_coverage_output_folder, final_output_folder, summary_output_folder,
    num_threads, qualified_quality_phred, reads_to_process,
    min_identity, min_coverage, reference_genome, skip_read_threshold
)

fastp_filtering.filter_reads()
