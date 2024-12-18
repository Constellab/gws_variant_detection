# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, ConfigSpecs, File, InputSpec, Folder, ShellProxy,
                      InputSpecs, IntParam, OutputSpec, OutputSpecs, Task, TableExporter,
                      TaskInputs, TaskOutputs, task_decorator, StrParam, Table , TableImporter , ResourceSet, TaskFileDownloader)


from ..base_env.Haplotypes_env import HaplotypesShellProxyHelper


@task_decorator("VariantsHaplotypesMildiou", human_name="VariantsHaplotypesMildiou",
                short_description="Haplotype Variant Calling")
class VariantsHaplotypesMildiou(Task):
    """
   This task enables systematic matrix construction of co-occurring variants.
    """

    input_specs = InputSpecs({
        'bam_folder':
        InputSpec(
            Folder, human_name="Folder containing the different alignment files"),
        'known_variants_file':
        InputSpec(
            Table, human_name="Table containing the result of filtered variant calling"),
        'ref_genome':
        InputSpec(
            File, human_name="File containing gene reference sequence")
    })

    output_specs = OutputSpecs({
        'counts_matrices_folder':
        OutputSpec(
            Folder, human_name="co-occurence matrix in term of counts",
            short_description="The frequency of reads of a target gene with multiple variants (mutation1 + mutations 2) in term of COUNTS"),

        'frequencies_matrices_folder':
        OutputSpec(
            Folder, human_name="co-occurence matrix in term of frequencies",
            short_description="The frequency of reads of a target gene with multiple variants (mutation1 + mutations 2) in term of FREQUENCY"),
        'resource_set':
        OutputSpec(
            ResourceSet, human_name="Count Co-occurence matrix",
            short_description="Multiple Variant Detection Tables in term of COUNTS"),
        'resource_set_FREQ':
        OutputSpec(
            ResourceSet, human_name="FREQ Co-occurence matrix",
            short_description="Multiple Variant Detection Tables in term of FREQUENCIES"),

        'wildtype_frequencies':OutputSpec(Table, human_name="wildtype_frequencies", short_description="Wildetype frequency per barcode"),

    })

    config_specs: ConfigSpecs = {
        "gene_name":
        StrParam(allowed_values=["Cytb","OSBP","Tubulin","ATPase"], default_value="Cytb",
                 short_description="Targeted gene")
    }

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_variant_haplotypes_mildiou.py"
    )

    shell_proxy: ShellProxy = None


    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        # retrive the input table
        bam_folder: Folder = inputs['bam_folder']
        known_variants_table: Table = inputs['known_variants_file']
        ref_genome: File = inputs['ref_genome']
        gene_name = params["gene_name"]

        # Setup the minipileup executable
        downloader = TaskFileDownloader(
            brick_name=VariantsHaplotypesMildiou.get_brick_name(),
            message_dispatcher=self.message_dispatcher)
        minipileup_folder = downloader.download_file_if_missing(
            url='https://storage.sbg.cloud.ovh.net/v1/AUTH_a0286631d7b24afba3f3cdebed2992aa/public/gws_variant_detection/minipileup-1.0.zip',
            filename='minipileup', decompress_file=True)

        minipileup_folder_v1 = os.path.join(os.path.join(minipileup_folder, 'minipileup-1.0'))

        #  path including the directory if exists
        minipileup_file_path = os.path.join(minipileup_folder_v1, 'minipileup')

        if not os.path.isfile(minipileup_file_path):
            # If the executable does not exist, navigate to the directory and run the make command
            shell_proxy = HaplotypesShellProxyHelper.create_proxy(self.message_dispatcher, minipileup_folder_v1)

            # Execute the make command
            make_command = "make"
            result = shell_proxy.run(make_command, shell_mode=True)
            if result:
                raise Exception("Error during build of minipileup, please check logs")

        # Execute the command
        shell_proxy = HaplotypesShellProxyHelper.create_proxy(self.message_dispatcher)
        output_folder_path = os.path.join(shell_proxy.working_dir, 'output_folder')
        os.makedirs(output_folder_path)


        counts_matrices_folder_path = os.path.join(shell_proxy.working_dir, 'counts_matrices_folder')
        os.makedirs(counts_matrices_folder_path)

        frequencies_matrices_folder_path = os.path.join(shell_proxy.working_dir, 'frequencies_matrices_folder')
        os.makedirs(frequencies_matrices_folder_path)


        file: File = TableExporter.call(known_variants_table, {
            'file_format': 'tsv',
            'delimiter' : 'tab',
            'write_header': True,
            'write_index' : False
        })


        # call python file
        cmd = f"python3 {self.python_file_path} {bam_folder.path} {ref_genome.path } {str(gene_name)} {minipileup_file_path} {output_folder_path} {counts_matrices_folder_path} {frequencies_matrices_folder_path} {file.path}"
        result = shell_proxy.run(cmd, shell_mode=True)

        if result != 0:
            raise Exception("Error during command, please check logs")
        self.shell_proxy = shell_proxy

        # display outputs :
        counts_matrices_folder = Folder(counts_matrices_folder_path)
        frequencies_matrices_folder = Folder(frequencies_matrices_folder_path)

        ## create a ressource_set
        resource_set = ResourceSet()
        for filename in os.listdir(counts_matrices_folder_path):
            file_path = os.path.join(counts_matrices_folder_path, filename)
            if os.path.isfile(file_path) and filename.endswith('.csv'):
                # Use TableImporter.call to read CSV files and add them to the resource set
                table = TableImporter.call(File(file_path))
                resource_set.add_resource(table, filename)

        ## create a ressource_set_freq
        resource_set_FREQ = ResourceSet()
        for filename in os.listdir(frequencies_matrices_folder_path):
            file_path = os.path.join(frequencies_matrices_folder_path, filename)
            if os.path.isfile(file_path) and filename.endswith('.csv'):
                # Use TableImporter.call to read CSV files and add them to the resource set
                table = TableImporter.call(File(file_path))
                resource_set_FREQ.add_resource(table, filename)


        ## Wildtype frequency per barcode :
        wildtype_frequencies_path = os.path.join(counts_matrices_folder_path, 'wildtype_frequencies.csv')
        wildtype_frequencies_table = TableImporter.call(File(wildtype_frequencies_path), {'delimiter': ',', 'header': 0, 'file_format': 'csv'})

        # Return the outputs
        return {
            'counts_matrices_folder': counts_matrices_folder,
            'frequencies_matrices_folder':frequencies_matrices_folder,
            'resource_set':resource_set,
            'resource_set_FREQ':resource_set_FREQ,
            'wildtype_frequencies':wildtype_frequencies_table
        }

    def run_after_task(self):
        if self.shell_proxy:
            self.shell_proxy.clean_working_dir()