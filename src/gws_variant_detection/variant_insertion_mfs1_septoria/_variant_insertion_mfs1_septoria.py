import os
import pysam
import subprocess
import pandas as pd
import re  # Import re for regex operations
import sys

class BAMClassifier:
    def __init__(self, bam_dir, output_dir, insertion_sizes_dir, reference_fasta,
                 type_1_tolerance=100, type2_338_tolerance=100, type2_369_tolerance=100,
                 type_3_tolerance=50, ref_tolerance=50):
        self.bam_dir = bam_dir
        self.output_dir = output_dir
        self.insertion_sizes_dir = insertion_sizes_dir
        self.reference_fasta = reference_fasta

        # Set specific tolerances for each classification type
        self.type_1_tolerance = type_1_tolerance
        self.type2_338_tolerance = type2_338_tolerance
        self.type2_369_tolerance = type2_369_tolerance
        self.type_3_tolerance = type_3_tolerance
        self.ref_tolerance = ref_tolerance

        # Set the adjusted reference size based on observed read length distribution
        self.ref_size = 552  # Mean read length for reference reads

    def classify_samples(self):
        # List BAM files in the directory
        bam_files = [f for f in os.listdir(self.bam_dir) if f.endswith(".bam")]
        summary_data = []

        # Create a folder for insertion sizes
        os.makedirs(self.insertion_sizes_dir, exist_ok=True)

        for bam_file in bam_files:
            bam_path = os.path.join(self.bam_dir, bam_file)
            sample_name = os.path.splitext(bam_file)[0]

            # Classify reads in the BAM file
            classified_reads, read_lengths = self._classify_reads(bam_path, sample_name)

            # Write classified reads to new BAM files
            for res_type, reads in classified_reads.items():
                bam_output_path = os.path.join(self.output_dir, f"{sample_name}_{res_type}.bam")
                self._write_bam(reads, bam_path, bam_output_path)

            # Calculate the total reads and frequencies for each classification type
            total_reads = sum(len(reads) for reads in classified_reads.values())
            summary_data.append({
                "Sample": sample_name,
                "Total_Reads": total_reads,
                "Type_1_Count": len(classified_reads["Type_1"]),
                "Type_1_Frequency": round((len(classified_reads["Type_1"]) / total_reads) * 100, 2) if total_reads > 0 else 0,
                "Type2_338_Count": len(classified_reads["Type2_338"]),
                "Type2_338_Frequency": round((len(classified_reads["Type2_338"]) / total_reads) * 100, 2) if total_reads > 0 else 0,
                "Type2_369_Count": len(classified_reads["Type2_369"]),
                "Type2_369_Frequency": round((len(classified_reads["Type2_369"]) / total_reads) * 100, 2) if total_reads > 0 else 0,
                "Type_3_Count": len(classified_reads["Type_3"]),
                "Type_3_Frequency": round((len(classified_reads["Type_3"]) / total_reads) * 100, 2) if total_reads > 0 else 0,
                "Ref_Similar_Count": len(classified_reads["Ref_Similar"]),
                "Ref_Similar_Frequency": round((len(classified_reads["Ref_Similar"]) / total_reads) * 100, 2) if total_reads > 0 else 0,
                "Other_Count": len(classified_reads["Other"]),
                "Other_Frequency": round((len(classified_reads["Other"]) / total_reads) * 100, 2) if total_reads > 0 else 0
            })

            # Save read lengths and classifications to a text file
            lengths_path = os.path.join(self.insertion_sizes_dir, f"{sample_name}_read_lengths.txt")
            with open(lengths_path, "w") as f:
                f.write("Sample\tRead_Name\tRead_Length\tClassification\n")
                for sample, read_name, length, classification in read_lengths:
                    f.write(f"{sample}\t{read_name}\t{length}\t{classification}\n")

        # Save summary data to CSV
        summary_df = pd.DataFrame(summary_data)

        # Function to extract the numeric barcode from the sample name
        def extract_barcode_number(sample_name):
            match = re.search(r'barcode(\d+)', sample_name)
            if match:
                return int(match.group(1))  # Extract the numeric barcode part as an integer
            return float('inf')  # If no barcode is found, put it at the end when sorting

        # Sort the DataFrame by the extracted barcode number
        summary_df['Barcode_Num'] = summary_df['Sample'].apply(extract_barcode_number)
        summary_df = summary_df.sort_values('Barcode_Num').drop('Barcode_Num', axis=1)

        # Write sorted DataFrame to a CSV file
        summary_csv_path = os.path.join(self.insertion_sizes_dir, "classification_summary.csv")
        summary_df.to_csv(summary_csv_path, index=False)

    def _classify_reads(self, bam_path, sample_name):
        classified_reads = {
            "Type_1": [],
            "Type2_338": [],
            "Type2_369": [],
            "Type_3": [],
            "Ref_Similar": [],
            "Other": []
        }
        read_lengths = []

        with pysam.AlignmentFile(bam_path, "rb") as bam_file:
            for read in bam_file:
                # Skip secondary or supplementary alignments
                if read.is_secondary or read.is_supplementary:
                    continue

                read_length = read.query_length
                res_type = self._classify_by_read_length(read_length)
                classified_reads[res_type].append(read)
                read_lengths.append((sample_name, read.query_name, read_length, res_type))

        return classified_reads, read_lengths

    def _classify_by_read_length(self, read_length):
        # Define expected sizes for each classification type based on reference size + insertion
        type_1_size = self.ref_size + 519
        type2_338_size = self.ref_size + 338
        type2_369_size = self.ref_size + 369
        type_3_size = self.ref_size + 150

        # Classification based on the calculated read length and specified tolerance for each type
        if (type_1_size - self.type_1_tolerance) <= read_length <= (type_1_size + self.type_1_tolerance):
            return "Type_1"
        elif (type2_338_size - self.type2_338_tolerance) <= read_length <= (type2_338_size + self.type2_338_tolerance):
            return "Type2_338"
        elif (type2_369_size - self.type2_369_tolerance) <= read_length <= (type2_369_size + self.type2_369_tolerance):
            return "Type2_369"
        elif (type_3_size - self.type_3_tolerance) <= read_length <= (type_3_size + self.type_3_tolerance):
            return "Type_3"
        elif (self.ref_size - self.ref_tolerance) <= read_length <= (self.ref_size + self.ref_tolerance):
            return "Ref_Similar"
        else:
            return "Other"

    def _write_bam(self, reads, original_bam_path, bam_output_path):
        # Write the classified reads to a new BAM file
        os.makedirs(os.path.dirname(bam_output_path), exist_ok=True)  # Ensure the directory exists
        with pysam.AlignmentFile(original_bam_path, "rb") as original_bam:
            with pysam.AlignmentFile(bam_output_path, "wb", header=original_bam.header) as out_bam:
                for read in reads:
                    out_bam.write(read)

        # Sort and index the BAM file
        sorted_bam_path = bam_output_path.replace(".bam", ".sorted.bam")
        subprocess.run(["samtools", "sort", "-o", sorted_bam_path, bam_output_path], check=True)
        subprocess.run(["samtools", "index", sorted_bam_path], check=True)

        # Remove the unsorted BAM file
        os.remove(bam_output_path)

# Example usage
bam_dir = sys.argv[1]
output_dir = sys.argv[2]
insertion_sizes_dir = sys.argv[3]
reference_fasta = sys.argv[4]  # Path to your reference genome

# Set individual tolerance values for each classification type
classifier = BAMClassifier(
    bam_dir, output_dir, insertion_sizes_dir, reference_fasta,
    type_1_tolerance=100, type2_338_tolerance=40, type2_369_tolerance=40,
    type_3_tolerance=100, ref_tolerance=50
)

classifier.classify_samples()
