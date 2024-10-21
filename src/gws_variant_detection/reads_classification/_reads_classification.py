import os
import gzip
import mappy as mp
import pysam
import sys

class ReadClassifier:
    def __init__(self, input_dir, output_dir, merged_reference_fasta, identity_threshold):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.merged_reference_fasta = merged_reference_fasta
        self.identity_threshold = identity_threshold

    def classify_reads_in_folder(self, folder_path):
        # Extract the name of the folder
        input_folder_name = os.path.basename(os.path.normpath(folder_path))

        # Load the merged reference sequence and create an aligner
        aligner = mp.Aligner(self.merged_reference_fasta, preset="map-ont")
        if not aligner:
            raise ValueError("Failed to load the merged reference sequence.")

        # Process each .fastq.gz file in the current folder
        for fastq_file in os.listdir(folder_path):
            if fastq_file.endswith(".fastq.gz"):
                with pysam.FastxFile(os.path.join(folder_path, fastq_file)) as infile:
                    for entry in infile:
                        for hit in aligner.map(entry.sequence):
                            # Calculate the identity as a percentage
                            identity = hit.mlen / hit.blen * 100
                            if identity >= self.identity_threshold:
                                gene_name = hit.ctg  # Directly use the reference name (e.g., tubulin)
                                # Define gene-specific output directory
                                gene_output_dir = os.path.join(self.output_dir, gene_name)
                                os.makedirs(gene_output_dir, exist_ok=True)  # Ensure the directory exists

                                # Define output file path within the gene-specific directory
                                output_file_path = os.path.join(gene_output_dir, f"{input_folder_name}-{gene_name}.fastq.gz")

                                # Open the file for writing, append mode
                                with gzip.open(output_file_path, 'at') as output_file:
                                    output_file.write(f"@{entry.name}\n{entry.sequence}\n+\n{entry.quality}\n")
                                break  # Exit loop after first successful classification

    def classify_reads(self):
        # Ensure the output directory exists
        os.makedirs(self.output_dir, exist_ok=True)

        # Process each subdirectory within the input directory
        for root, dirs, files in os.walk(self.input_dir):
            for dir in dirs:
                self.classify_reads_in_folder(os.path.join(root, dir))

    # Define input and output directories, and the path to the merged reference sequence
input_dir = sys.argv[1]
output_dir = sys.argv[3]
merged_reference_fasta = sys.argv[2]
identity_threshold = 60  # Example threshold

    # Initialize and run the classifier
classifier = ReadClassifier(input_dir, output_dir, merged_reference_fasta, identity_threshold)
classifier.classify_reads()
