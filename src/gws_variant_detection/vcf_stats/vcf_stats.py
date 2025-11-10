# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, ConfigSpecs, InputSpec, Folder,File,
                      InputSpecs, OutputSpec, OutputSpecs, Task,
                      TaskInputs, TaskOutputs, task_decorator, Table , TableImporter)


from ..base_env.Variants_env import VariantsShellProxyHelper


@task_decorator("VCFstats", human_name="VariantCallingStats",
                short_description="Statistical variant calling file")
class VariantStats(Task):
    """
    This task automates the generation and parsing of statistics from Variant Call Format (VCF) files using the bcftools stats command.
    It provides an overview of various types of genetic variants detected in the input VCF files, including :
    SNPs (single nucleotide polymorphisms) , indels (insertions/deletions) , MNP (multi-nucleotide polymorphisms e.g. a dinucleotide substitution ) ,
    and complex variants (combination of SNPs and indels at a single position).

    """

    input_specs = InputSpecs({
        'input_folder':
        InputSpec(
            Folder, human_name="Folder containing the different normalized variant calling files"),
    })

    output_specs = OutputSpecs({
        'Statistical_vcf_table': OutputSpec(Table, human_name="VCF stats", short_description="An overview of the numbers of the different detected variants.")
    })

    config_specs = ConfigSpecs({})

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_vcf_stats.py"
    )

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        # retrive the input table
        input_folder: Folder = inputs['input_folder']

        # Execute the command
        shell_proxy = VariantsShellProxyHelper.create_proxy(self.message_dispatcher)
        output_path = os.path.join(shell_proxy.working_dir, 'vcf_stat')
        os.makedirs(output_path)

        # call python file
        cmd = f"python3 {self.python_file_path} {input_folder.path} {output_path}"
        shell_proxy.run(cmd, shell_mode=True)

        # display outputs :

        # Read the known and unknown variants from the output path
        statistical_variants_file_path = os.path.join(output_path, 'summary_stats.tsv')

        statistical_variants_table = TableImporter.call(File(statistical_variants_file_path), {'delimiter': 'tab', 'header': 0, 'file_format': 'tsv'})

        # Return the outputs
        return {
            'Statistical_vcf_table': statistical_variants_table
        }