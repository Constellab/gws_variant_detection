# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, ConfigSpecs, File, InputSpec, Folder,ShellProxy,
                      InputSpecs, IntParam, OutputSpec, OutputSpecs, Task,
                      TaskInputs, TaskOutputs, task_decorator, StrParam, Table , TableImporter)


from ..base_env.Variants_env import VariantsShellProxyHelper


@task_decorator("SeptoriaInsertionsMfs1", human_name="SeptoriaInsertionsMfs1",
                short_description="Determination of insertion size of the MFS1 gene in Septoria.")
class SeptoriaInsertionsMfs1(Task):
    """
    This script classifies reads from BAM files based on their insertion sizes relative to a reference genome and generates categorized output BAM files and summary statistics. It calculates insertion sizes, classifies reads into defined types, and outputs the results in CSV and BAM formats
    """

    input_specs = InputSpecs({
        'input_folder':
        InputSpec(
            Folder, human_name="Folder containing the different alignment files"),
        'input_file':
        InputSpec(
            File, human_name="File containing MFS1 gene reference sequence")
    })

    output_specs = OutputSpecs({
        'output_bam_path':
        OutputSpec(
            Folder, human_name="Per insertions alignment files",
            short_description="A folder containing the different alignment files based on insertion size"),
        'output_insertion_size_path':
        OutputSpec(
            Folder, human_name="Per read insertion size in each sample",
            short_description="A folder containing the different insertion size in each read"),
        'classification_summary': OutputSpec(Table, human_name="Classification summary", short_description="Table contains a summary of the read classifications for each sample"),

    })

    config_specs: ConfigSpecs = {
        "mdr_tolerence": IntParam(default_value=20, min_value=0, short_description="This parameter defines the acceptable variation (in base pairs) around the expected insertion sizes for classifying reads. "),
        "ref_tolerance": IntParam(default_value=20, min_value=0, short_description="This parameter specifies the tolerance (in base pairs) for classifying reads as Ref_Similar indicating how much deviation from the reference size is allowed for a read to be considered similar to the reference.")

    }

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_variant_insertion_mfs1.py"
    )

    shell_proxy: ShellProxy = None

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        # retrive the input table
        input_folder: Folder = inputs['input_folder']
        input_file: File = inputs['input_file']
        mdr_tolerence = params["mdr_tolerence"]
        ref_tolerance = params["ref_tolerance"]

        # Execute the command
        shell_proxy = VariantsShellProxyHelper.create_proxy(self.message_dispatcher)
        alignment_result_path = os.path.join(shell_proxy.working_dir, 'Alignment_bam_files')
        os.makedirs(alignment_result_path)

        per_read_insertion_size_result_path = os.path.join(shell_proxy.working_dir, 'per_read_insertion_size')
        os.makedirs(per_read_insertion_size_result_path)


        # call python file
        cmd = f"python3 {self.python_file_path} {input_folder.path} {alignment_result_path} {per_read_insertion_size_result_path} {input_file.path} {int(mdr_tolerence)} {int(ref_tolerance)}"
        shell_proxy.run(cmd, shell_mode=True)
        self.shell_proxy = shell_proxy

        # display outputs :
        # 1_normalized vcf
        alignment_bam_files = Folder(alignment_result_path)
        per_read_insertion_size_result_files = Folder(per_read_insertion_size_result_path)

        # Read the known and unknown variants from the output path
        classification_summary_file_path = os.path.join(per_read_insertion_size_result_path, 'classification_summary.csv')

        classification_summary_table = TableImporter.call(File(classification_summary_file_path), {'delimiter': ',', 'header': 0, 'file_format': 'csv'})

        # Return the outputs
        return {
            'output_bam_path': alignment_bam_files,
            'output_insertion_size_path': per_read_insertion_size_result_files,
            'classification_summary': classification_summary_table,

        }

    def run_after_task(self):
        if self.shell_proxy:
            self.shell_proxy.clean_working_dir()