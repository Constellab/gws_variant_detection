import sys
import os
import subprocess
import gzip
import random
import pandas as pd


class FastpFiltering:

    def __init__(
        self,
        input_folder: str,
        fastp_output_folder: str,
        final_output_folder: str,
        summary_output_folder: str,
        num_threads: int,
        qualified_quality_phred: int,
        reads_to_process: int,
        skip_read_threshold: int,
        min_read_length: int,
    ) -> None:
        self.input_folder = input_folder
        self.fastp_output_folder = fastp_output_folder
        self.final_output_folder = final_output_folder
        self.summary_output_folder = summary_output_folder
        self.num_threads = num_threads
        self.qualified_quality_phred = qualified_quality_phred
        self.reads_to_process = reads_to_process
        self.skip_read_threshold = skip_read_threshold
        self.min_read_length = min_read_length

        os.makedirs(self.fastp_output_folder, exist_ok=True)
        os.makedirs(self.final_output_folder, exist_ok=True)
        os.makedirs(self.summary_output_folder, exist_ok=True)


    def filter_reads(self) -> None:
        """Process every ``*.fastq.gz`` file in *input_folder*."""
        summary_records: list[dict] = []

        fastq_files = [f for f in os.listdir(self.input_folder) if f.endswith(".fastq.gz")]
        if not fastq_files:
            print("[Error] No .fastq.gz files found in input folder", file=sys.stderr)
            return

        for fq in fastq_files:
            src = os.path.join(self.input_folder, fq)
            fastp_out = os.path.join(self.fastp_output_folder, fq)
            final_out = os.path.join(self.final_output_folder, fq)

            # 1) filtering via *chopper*
            cmd = (
                f"gunzip -c {src} | "
                f"chopper -q {self.qualified_quality_phred} --threads {self.num_threads} | "
                f"gzip > {fastp_out}"
            )
            if subprocess.run(cmd, shell=True).returncode != 0:
                print(f"[Error] chopper failed on {fq}", file=sys.stderr)
                continue

            # 2) Length filter
            valid_reads: list[tuple[str, str, str]] = []  # (name, seq, qual)
            with gzip.open(fastp_out, "rt") as handle:
                lines = handle.readlines()

            for i in range(0, len(lines), 4):
                name = lines[i].rstrip()
                seq = lines[i + 1].rstrip()
                qual = lines[i + 3].rstrip()

                if len(seq) >= self.min_read_length:
                    valid_reads.append((name, seq, qual))

            # 3) Random subsampling
            if len(valid_reads) > self.reads_to_process:
                rng = random.Random(123)
                final_reads = rng.sample(valid_reads, self.reads_to_process)
            else:
                final_reads = valid_reads

            # 4) Write final FASTQ
            with gzip.open(final_out, "wt") as fh:
                for name, seq, qual in final_reads:
                    fh.write(f"{name}\n{seq}\n+\n{qual}\n")

            # 5) Summary stats
            summary_records.append(
                {
                    "file": fq,
                    "reads_quality_length": len(valid_reads),
                    "reads_final": len(final_reads),
                    "phred_threshold": self.qualified_quality_phred,
                    "min_read_length": self.min_read_length,
                    "reads_target": self.reads_to_process,
                }
            )

        # 6) Global CSV
        if summary_records:
            df = pd.DataFrame(summary_records)
            csv_path = os.path.join(self.summary_output_folder, "filtering_summary.csv")
            df.to_csv(csv_path, index=False)
            print(f"[Info] Summary written to {csv_path}")
        else:
            print("[Warning] No summary generated – nothing processed.")


if __name__ == "__main__":
    USAGE = (
        "Usage: python script.py "
        "<input_folder> <fastp_output_folder> <final_output_folder> <summary_output_folder> "
        "<num_threads> <phred_threshold> <reads_to_process>"
        "<skip_read_threshold> <min_read_length>"
    )

    (
        input_folder,
        fastp_output_folder,
        final_output_folder,
        summary_output_folder,
        num_threads,
        phred_threshold,
        reads_to_process,
        skip_read_threshold,
        min_read_length,
    ) = sys.argv[1:10]

    # Type conversions
    num_threads = int(num_threads)
    phred_threshold = int(phred_threshold)
    reads_to_process = int(reads_to_process)
    skip_read_threshold = int(skip_read_threshold)
    min_read_length = int(min_read_length)

    pipeline = FastpFiltering(
        input_folder=input_folder,
        fastp_output_folder=fastp_output_folder,
        final_output_folder=final_output_folder,
        summary_output_folder=summary_output_folder,
        num_threads=num_threads,
        qualified_quality_phred=phred_threshold,
        reads_to_process=reads_to_process,
        skip_read_threshold=skip_read_threshold,
        min_read_length=min_read_length,
    )

    pipeline.filter_reads()
