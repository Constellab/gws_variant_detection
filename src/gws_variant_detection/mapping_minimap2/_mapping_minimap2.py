import os
import subprocess
import sys

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
        subprocess.run(["minimap2", "-d", f"{self.ref_genome}.mmi", self.ref_genome], check=True)
        subprocess.run(["samtools", "faidx", self.ref_genome], check=True)

    def map_reads(self):
        file_extension = ".fastq.gz"
        platform = "Nanopore"
        read_files = [f for f in os.listdir(self.fastq_dir) if f.endswith(file_extension)]

        RG_COUNTER = 0
        for read_file in read_files:
            bname = os.path.splitext(read_file)[0]

            RG_COUNTER += 1
            RGID = f"RG{RG_COUNTER}"
            RGSN = bname
            readgroup = f"@RG\\tID:{RGID}\\tPL:{platform}\\tPU:unit1\\tLB:lib1\\tSM:{RGSN}"

            print(f"Mapping read file: {read_file}")

            output_bam = os.path.join(self.output_dir, f"{bname}.bam")

            # Mapping avec minimap2 et conversion en BAM avec samtools
            minimap_cmd = ["minimap2", "-t", str(self.threads), "-R", readgroup, "-ax", "map-ont", self.ref_genome, os.path.join(self.fastq_dir, read_file)]
            samtools_cmd = ["samtools", "view", "-bS", "-"]

            minimap_process = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE)
            samtools_process = subprocess.Popen(samtools_cmd, stdin=minimap_process.stdout, stdout=subprocess.PIPE)
            with open(output_bam, 'wb') as f:
                for data in samtools_process.stdout:
                    f.write(data)
            minimap_process.stdout.close()
            samtools_process.communicate()

            # Filtrer pour ne conserver que les alignements primaires
            # -F 2304 permet d'exclure les alignements secondaires (0x100) et complémentaires (0x800)
            primary_bam = os.path.join(self.output_dir, f"{bname}.primary.bam")
            subprocess.run(["samtools", "view", "-b", "-F", "2304", "-o", primary_bam, output_bam], check=True)

            sorted_bam = os.path.join(self.output_dir, f"{bname}.sorted.bam")
            subprocess.run(["samtools", "sort", "-o", sorted_bam, primary_bam], check=True)
            subprocess.run(["samtools", "index", sorted_bam], check=True)

            os.remove(output_bam)
            os.remove(primary_bam)

    def preprocess_bams(self):
        for f in os.listdir(self.output_dir):
            if f.endswith(".bam"):
                new_name = f.replace(".fastq.sorted.bam", ".sorted.bam")
                os.rename(os.path.join(self.output_dir, f), os.path.join(self.output_dir, new_name))
            if f.endswith(".bai"):
                new_name = f.replace(".fastq.sorted.bam.bai", ".sorted.bam.bai")
                os.rename(os.path.join(self.output_dir, f), os.path.join(self.output_dir, new_name))

        print("Files have been renamed to match the desired format.")

fastq_dir = sys.argv[1]
ref_genome = sys.argv[2]
threads = int(sys.argv[3])
output_dir = sys.argv[4]

mapper = ReadsMappingMinimap(ref_genome, fastq_dir, threads, output_dir)
mapper.create_output_directory()
mapper.index_reference()
mapper.map_reads()
mapper.preprocess_bams()