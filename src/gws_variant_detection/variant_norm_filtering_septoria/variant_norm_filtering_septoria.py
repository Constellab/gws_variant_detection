# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, ConfigSpecs, File, InputSpec, Folder,ShellProxy,
                      InputSpecs, IntParam, OutputSpec, OutputSpecs, Task,
                      TaskInputs, TaskOutputs, task_decorator, StrParam, Table , TableImporter)


from ..base_env.Variants_env import VariantsShellProxyHelper


@task_decorator("VariantsFilteringSeptoria", human_name="VariantsFilteringSeptoria",
                short_description="Filtering the different detected variant belonging to SEPTORIA")
class VariantsFilteringSeptoria(Task):
    """
    This task is designed to normalize and process Variant Call Format (VCF) files obtained from next-generation sequencing experiments.
    It takes input VCF files containing genetic variant data, normalizes them, extracts relevant information, merges the data from multiple files, identifies known variants for specific genes, and separates them into known and unknown variants.
    """

    input_specs = InputSpecs({
        'input_folder':
        InputSpec(
            Folder, human_name="Folder containing the different variant calling files"),
        'input_file':
        InputSpec(
            File, human_name="File containing gene reference sequence")
    })

    output_specs = OutputSpecs({
        'output_path':
        OutputSpec(
            Folder, human_name="Normalized variant calling files",
            short_description="A folder containing the different normalized variant calling files"),
        'known_variants': OutputSpec(Table, human_name="Known variants", short_description="A list of variants that have been previously identified and reported"),
        'unknown_variants': OutputSpec(Table, human_name="Unknow variants", short_description="A list of variant of uncertain (or unknown)"),

    })

    config_specs: ConfigSpecs = {
        "gene_name":
        StrParam(allowed_values=["CYP51","CytB","SdhB","SdhC","SdhD"], default_value="CytB",
                 short_description="Targeted gene"),
        "reads_to_process": IntParam(default_value=5000, min_value=3000, short_description="Number of reads to process")
    }

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_variant_norm_filtering_septoria.py"
    )

    shell_proxy: ShellProxy = None

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        # retrive the input table
        input_folder: Folder = inputs['input_folder']
        input_file: File = inputs['input_file']
        gene_name = params["gene_name"]
        reads_to_process = params["reads_to_process"]

        # Execute the command
        shell_proxy = VariantsShellProxyHelper.create_proxy(self.message_dispatcher)
        normalized_result_path = os.path.join(shell_proxy.working_dir, 'Normalized_vcf_files')
        os.makedirs(normalized_result_path)

        known_unknown_result_path = os.path.join(shell_proxy.working_dir, 'known_unknown_files')
        os.makedirs(known_unknown_result_path)


        # call python file
        cmd = f"python3 {self.python_file_path} {input_folder.path} {input_file.path } {str(gene_name)} {int(reads_to_process)} {normalized_result_path} {known_unknown_result_path}"
        shell_proxy.run(cmd, shell_mode=True)
        self.shell_proxy = shell_proxy

        # display outputs :
        # 1_normalized vcf
        normalized_vcf_files = Folder(normalized_result_path)

        # Read the known and unknown variants from the output path
        known_variants_file_path = os.path.join(known_unknown_result_path, 'known_variants_filtered.tsv')
        unknown_variants_file_path = os.path.join(known_unknown_result_path, 'unknown_variants.tsv')

        known_variants_table = TableImporter.call(File(known_variants_file_path), {'delimiter': 'tab', 'header': 0, 'file_format': 'tsv'})
        unknown_variants_table = TableImporter.call(File(unknown_variants_file_path), {'delimiter': 'tab', 'header': 0, 'file_format': 'tsv'})

        # Return the outputs
        return {
            'output_path': normalized_vcf_files,
            'known_variants': known_variants_table,
            'unknown_variants': unknown_variants_table
        }

    def run_after_task(self):
        if self.shell_proxy:
            self.shell_proxy.clean_working_dir()