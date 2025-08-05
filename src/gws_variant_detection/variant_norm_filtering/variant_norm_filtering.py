# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, ConfigSpecs, File, InputSpec, Folder, ShellProxy,
                      InputSpecs, IntParam, OutputSpec, OutputSpecs, Task,
                      TaskInputs, TaskOutputs, task_decorator, StrParam, Table, TableImporter, ResourceSet, TaskFileDownloader)


from ..base_env.Variants_env import VariantsShellProxyHelper


@task_decorator("VariantsFiltering", human_name="VariantsFiltering",
                short_description="Filtering the different detected variant")
class VariantsFiltering(Task):
    """
    This task is designed to normalize and process Variant Call Format (VCF) files obtained from next-generation sequencing experiments.
    It takes input VCF files containing genetic variant data, normalizes them, extracts relevant information, merges the data from multiple files, identifies known variants for specific genes, and separates them into known and unknown variants.
    """

    input_specs = InputSpecs({
        'input_folder':
        InputSpec(
            Folder, human_name="Folder containing the different BAM files"),
        'input_file':
        InputSpec(
            File, human_name="File containing gene reference sequence"),
    })

    output_specs = OutputSpecs({
        'read_variants_table':OutputSpec(Table, human_name="read variants table", short_description="A table containing the known variants detected per read and per barcode"),
        'filtered_reads_count':OutputSpec(Table, human_name="filtered reads count", short_description="A table containing the number of reads per barcode after alignment"),
        'known_variants_filtered': OutputSpec(Table, human_name="Known variants", short_description="A list of variants that have been previously identified and reported"),
        'unknown_variants': OutputSpec(Table, human_name="Unknow variants", short_description="A list of variant of uncertain (or unknown)"),
        'wildtype_frequencies':OutputSpec(Table, human_name="wildtype_frequencies", short_description="Wildetype frequency per barcode"),

    })

    config_specs: ConfigSpecs = {
        "gene_name":
        StrParam(allowed_values=["Cytb", "OSBP", "Tubulin", "ATPase"], default_value="Cytb",
                 short_description="Targeted gene"),
        "reads_to_process": IntParam(default_value=5000, min_value=3000, short_description="Number of reads to process"),
        "num_threads": IntParam(default_value=30, min_value=2, short_description="Number of threads "),

    }

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_variant_norm_filtering.py"
    )

    shell_proxy: ShellProxy = None

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        # retrive the input table
        input_folder: Folder = inputs['input_folder']
        input_file: File = inputs['input_file']
        gene_name = params["gene_name"]
        reads_to_process = params["reads_to_process"]
        num_threads = params["num_threads"]

        # Setup the minipileup executable
        downloader = TaskFileDownloader(
            brick_name=VariantsFiltering.get_brick_name(),
            message_dispatcher=self.message_dispatcher)
        minipileup_folder = downloader.download_file_if_missing(
            url='https://storage.sbg.cloud.ovh.net/v1/AUTH_a0286631d7b24afba3f3cdebed2992aa/public/gws_variant_detection/minipileup-1.0.zip',
            filename='minipileup', decompress_file=True)

        minipileup_folder_v1 = os.path.join(os.path.join(minipileup_folder, 'minipileup-1.0'))

        #  path including the directory if exists
        minipileup_file_path = os.path.join(minipileup_folder_v1, 'minipileup')

        if not os.path.isfile(minipileup_file_path):
            # If the executable does not exist, navigate to the directory and run the make command
            shell_proxy = VariantsShellProxyHelper.create_proxy(self.message_dispatcher, minipileup_folder_v1)

            # Execute the make command
            make_command = "make"
            result = shell_proxy.run(make_command, shell_mode=True)
            if result:
                raise Exception("Error during build of minipileup, please check logs")


        # Execute the command
        shell_proxy = VariantsShellProxyHelper.create_proxy()
        known_unknown_result_path = os.path.join(shell_proxy.working_dir, 'known_unknown_files')
        os.makedirs(known_unknown_result_path)


        # call python file
        cmd = f"python3 {self.python_file_path} {input_folder.path} {input_file.path } {str(gene_name)} {minipileup_file_path} {known_unknown_result_path} {int(num_threads)} {int(reads_to_process)}"
        shell_proxy.run(cmd, shell_mode=True)
        self.shell_proxy = shell_proxy

        # display outputs :
        read_variants_table_file_path = os.path.join(known_unknown_result_path, 'read_variants_table.csv')
        known_variants_file_path = os.path.join(known_unknown_result_path, 'known_variants_filtered.tsv')
        unknown_variants_file_path = os.path.join(known_unknown_result_path, 'unknown_variants.tsv')
        filtered_reads_count_file_path = os.path.join(known_unknown_result_path, 'filtered_reads_count.csv')
        wildtype_frequencies_file_path = os.path.join(known_unknown_result_path, 'wildtype_frequencies.csv')

        # Read the known and unknown variants from the output path

        read_variants_table_result = TableImporter.call(
            File(read_variants_table_file_path),
            {'delimiter': ',', 'header': 0, 'file_format': 'csv'})
        known_variants_table = TableImporter.call(
            File(known_variants_file_path),
            {'delimiter': 'tab', 'header': 0, 'file_format': 'tsv'})
        unknown_variants_table = TableImporter.call(
            File(unknown_variants_file_path),
            {'delimiter': 'tab', 'header': 0, 'file_format': 'tsv'})
        filtered_reads_count_table = TableImporter.call(
            File(filtered_reads_count_file_path),
            {'delimiter': ',', 'header': 0, 'file_format': 'csv'})
        wildtype_frequencies_table = TableImporter.call(
            File(wildtype_frequencies_file_path),
            {'delimiter': ',', 'header': 0, 'file_format': 'csv'})
        # Return the outputs
        return {
            'read_variants_table': read_variants_table_result,
            'known_variants_filtered': known_variants_table,
            'unknown_variants': unknown_variants_table,
            'filtered_reads_count': filtered_reads_count_table,
            'wildtype_frequencies':wildtype_frequencies_table
        }

    def run_after_task(self):
        if self.shell_proxy:
            self.shell_proxy.clean_working_dir()
