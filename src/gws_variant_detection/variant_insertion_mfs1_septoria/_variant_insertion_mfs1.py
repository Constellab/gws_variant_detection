import os
import pysam
import subprocess
import pandas as pd
import re  # Import re for regex operations
import sys

class BAMClassifier:
    def __init__(self, bam_dir, output_dir, insertion_sizes_dir, reference_fasta, tolerance=30, ref_tolerance=30):
        self.bam_dir = bam_dir
        self.output_dir = output_dir
        self.insertion_sizes_dir = insertion_sizes_dir
        self.tolerance = tolerance
        self.ref_tolerance = ref_tolerance
        self.reference_fasta = reference_fasta
        self.ref_size = self._get_reference_size()

    def _get_reference_size(self):
        # Calculate reference size using pysam
        with pysam.FastaFile(self.reference_fasta) as fasta:
            ref_length = 0
            for seq in fasta.references:
                ref_length = fasta.get_reference_length(seq)
                break
        print(f"Reference sequence size: {ref_length} bp")

        return ref_length

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
            classified_reads, insertion_sizes = self._classify_reads(bam_path, sample_name)

            # Write classified reads to new BAM files
            for res_type, reads in classified_reads.items():
                bam_output_path = os.path.join(self.output_dir, f"{sample_name}_{res_type}.bam")
                self._write_bam(reads, bam_path, bam_output_path)

            # Collect summary data
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

            # Save insertion sizes to the newly created insertion_sizes directory
            insertion_sizes_path = os.path.join(self.insertion_sizes_dir, f"{sample_name}_insertion_sizes.txt")
            with open(insertion_sizes_path, "w") as f:
                f.write("Sample\tRead_Name\tInsertion_Size\tClassification\n")
                for sample, read_name, size, classification in insertion_sizes:
                    f.write(f"{sample}\t{read_name}\t{size}\t{classification}\n")

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
        insertion_sizes = []

        with pysam.AlignmentFile(bam_path, "rb") as bam_file:
            for read in bam_file:
                # Skip secondary or supplementary alignments
                if read.is_secondary or read.is_supplementary:
                    continue

                insertion_size = self._calculate_insertion_size(read)
                res_type = self._classify_insertion(insertion_size)
                classified_reads[res_type].append(read)
                insertion_sizes.append((sample_name, read.query_name, insertion_size, res_type))

        return classified_reads, insertion_sizes

    def _calculate_insertion_size(self, read):
        # Calculate insertion size based on CIGAR string
        insertion_size = 0
        for (operation, length) in read.cigartuples:
            if operation == 1:  # '1' is the BAM_CINS code for insertion
                insertion_size += length
        return insertion_size

    def _classify_insertion(self, insertion_size):
        # Total size of the read considering the reference size
        total_size = self.ref_size + insertion_size

        # Classification based on total size
        if total_size >= (self.ref_size + 519) - self.tolerance and total_size <= (self.ref_size + 519) + self.tolerance:  # Type_1
            return "Type_1"
        elif total_size >= (self.ref_size + 338) - self.tolerance and total_size <= (self.ref_size + 338) + self.tolerance:  # Type2_338
            return "Type2_338"
        elif total_size >= (self.ref_size + 369) - self.tolerance and total_size <= (self.ref_size + 369) + self.tolerance:  # Type2_369
            return "Type2_369"
        elif total_size >= (self.ref_size + 150) - self.tolerance and total_size <= (self.ref_size + 150) + self.tolerance:  # Type_3
            return "Type_3"
        elif (self.ref_size - self.ref_tolerance) <= total_size <= (self.ref_size + self.ref_tolerance):  # Ref_Similar
            return "Ref_Similar"
        else:
            return "Other"

    def _write_bam(self, reads, original_bam_path, bam_output_path):
        # Write the classified reads to a new BAM file
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
tolerance = int(sys.argv[5])
ref_tolerance = int(sys.argv[6])

# Create the directories for output and insertion sizes
os.makedirs(output_dir, exist_ok=True)
os.makedirs(insertion_sizes_dir, exist_ok=True)

classifier = BAMClassifier(bam_dir, output_dir, insertion_sizes_dir, reference_fasta, tolerance=tolerance, ref_tolerance=ref_tolerance)
classifier.classify_samples()
