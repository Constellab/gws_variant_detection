import sys
import os
import subprocess

class ReadsMappingMinimap:
    def __init__(self, ref_genome, fastq_dir, threads, output_dir):
        self.ref_genome = ref_genome
        self.fastq_dir = fastq_dir
        self.threads = threads
        self.output_dir = output_dir

    def create_output_directory(self):
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    def index_reference(self):
        subprocess.run(["minimap2", "-d", f"{self.ref_genome}.mmi", self.ref_genome])
        subprocess.run(["samtools", "faidx", self.ref_genome])

    def map_reads(self):
        file_extension = ".fastq.gz"
        platform = "Nanopore"
        read_files = [f for f in os.listdir(self.fastq_dir) if f.endswith(file_extension)]

        RG_COUNTER = 0
        for read_file in read_files:
            bname = os.path.splitext(read_file)[0]

            RG_COUNTER += 1
            RGID = str(RG_COUNTER)
            RGSN = bname
            readgroup = f"@RG\\tID:{RGID}\\tPL:{platform}\\tPU:None\\tLB:None\\tSM:{RGSN}"

            print(f"read file: {read_file}")

            output_bam = os.path.join(self.output_dir, f"{bname}.raw.bam")
            minimap_cmd = ["minimap2", "-t", str(self.threads), "-R", readgroup, "-ax", "map-ont", self.ref_genome, os.path.join(self.fastq_dir, read_file)]
            samtools_cmd = ["samtools", "view", "-h", "-b", "-o", output_bam]

            minimap_process = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE)
            samtools_process = subprocess.Popen(samtools_cmd, stdin=minimap_process.stdout)
            minimap_process.stdout.close()
            samtools_process.communicate()

            # Rename output files to match the desired format
            os.rename(output_bam, os.path.join(self.output_dir, f"{bname}.raw.bam"))

    def preprocess_bams(self):
        allbams = [f for f in os.listdir(self.output_dir) if f.endswith(".raw.bam")]

        for raw_bam in allbams:
            bname = os.path.splitext(raw_bam)[0]

            filtered_bam = os.path.join(self.output_dir, f"{bname}.filtered.bam")
            sorted_bam_n = os.path.join(self.output_dir, f"{bname}.sorted.n.bam")
            fixmate_bam = os.path.join(self.output_dir, f"{bname}.fixmate.bam")
            sorted_bam_p = os.path.join(self.output_dir, f"{bname}.sorted.p.bam")
            dedup_bam = os.path.join(self.output_dir, f"{bname}.dedup.bam")
            dedup_bam_new_name = os.path.join(self.output_dir, f"{bname}.dedup.bam")

            samtools_view_cmd = ["samtools", "view", "-b", "-F", "0xc", os.path.join(self.output_dir, raw_bam), "-o", filtered_bam]
            samtools_sort_n_cmd = ["samtools", "sort", "-n", "-@", str(self.threads), filtered_bam, "-o", sorted_bam_n]
            samtools_fixmate_cmd = ["samtools", "fixmate", "-m", sorted_bam_n, fixmate_bam]
            samtools_sort_p_cmd = ["samtools", "sort", "-@", str(self.threads), fixmate_bam, "-o", sorted_bam_p]
            samtools_markdup_cmd = ["samtools", "markdup", "-r", "-@", str(self.threads), sorted_bam_p, dedup_bam_new_name]
            samtools_index_cmd = ["samtools", "index", dedup_bam_new_name]

            subprocess.run(samtools_view_cmd)
            subprocess.run(samtools_sort_n_cmd)
            subprocess.run(samtools_fixmate_cmd)
            subprocess.run(samtools_sort_p_cmd)
            subprocess.run(samtools_markdup_cmd)
            subprocess.run(samtools_index_cmd)

            os.remove(filtered_bam)
            os.remove(sorted_bam_n)
            os.remove(fixmate_bam)
            os.remove(sorted_bam_p)

            print("duplicates have been removed")

        # Rename files to match desired format
        for f in os.listdir(self.output_dir):
            if f.endswith(".raw.bam"):
                new_name = f.replace(".fastq.raw.bam", ".raw.bam")
                os.rename(os.path.join(self.output_dir, f), os.path.join(self.output_dir, new_name))
            elif f.endswith(".raw.dedup.bam"):
                new_name = f.replace(".fastq.raw.dedup.bam", ".dedup.bam")
                os.rename(os.path.join(self.output_dir, f), os.path.join(self.output_dir, new_name))
            elif f.endswith(".raw.dedup.bam.bai"):
                new_name = f.replace(".fastq.raw.dedup.bam.bai", ".dedup.bam.bai")
                os.rename(os.path.join(self.output_dir, f), os.path.join(self.output_dir, new_name))


fastq_dir = sys.argv[1]
ref_genome = sys.argv[2]
threads = int(sys.argv[3])
output_dir = sys.argv[4]
mapper = ReadsMappingMinimap(ref_genome, fastq_dir, threads, output_dir)
mapper.create_output_directory()
mapper.index_reference()
mapper.map_reads()
mapper.preprocess_bams()
