# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, ConfigSpecs, File, InputSpec, Folder, ShellProxy,
                      InputSpecs, IntParam, OutputSpec, OutputSpecs, Task,
                      TaskInputs, TaskOutputs, task_decorator,IntParam, Table,TableImporter)


from ..base_env.Filtering_env import FilteringShellProxyHelper

@task_decorator("ReadsQualityFilteringMfs1", human_name="NanoReadsFilteringMfs1",
                short_description="Elimination of low-quality reads and fixing the number of reads per sample.")
class ReadsQualityFilteringMfs1(Task):
    """
    This task automates the process of filtering and randomly sampling reads from sequencing data using chopper and seqtk respectively.
    The number of threads, a quality threshold, and the number of reads to process must be specified.
    """

    input_specs = InputSpecs({
        'input_folder':
        InputSpec(
            Folder, human_name="Folder containing fastq files"),
        'reference_fasta':
        InputSpec(
            File, human_name="Combined reference sequences in FASTA format")
    })

    output_specs = OutputSpecs({
        'final_filtered_reads':
        OutputSpec(
            Folder, human_name="Subset of filtered reads",
            short_description="This folder contains a subset of randomly filtered reads"),
        'number_reads_filtered':
        OutputSpec(
            Table, human_name="number of reads selected per barcode",
            short_description="This table contains the number of reads selected per barcode after filtration")
    })

    config_specs: ConfigSpecs = {
        "num_threads": IntParam(default_value=15, min_value=2, short_description="Number of threads "),
        "qualified_quality_phred": IntParam(default_value=10, min_value=4, short_description="qualified quality phred "),
        "reads_to_process": IntParam(default_value=5000, min_value=1000, short_description="number of reads to process"),
        "min_identity": IntParam(default_value=60, min_value=50, visibility=IntParam.PROTECTED_VISIBILITY, short_description=" Identity threshold"),
        "min_coverage": IntParam(default_value=40, min_value=30, visibility=IntParam.PROTECTED_VISIBILITY, short_description="Coverage threshold"),
        "skip_read_threshold": IntParam(default_value=100, min_value=1, human_name="skip_read_threshold", visibility=IntParam.PROTECTED_VISIBILITY,
            short_description="Threshold for skipping files with slightly fewer reads"),
        "min_read_length": IntParam(default_value=400, min_value=390, human_name="min_read_length", visibility=IntParam.PROTECTED_VISIBILITY,
            short_description="minimim read length")

    }

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_filtering_chopper_mfs1.py"
    )

    shell_proxy: ShellProxy = None

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        # retrive the input table
        input_folder: Folder = inputs['input_folder']
        reference_fasta: File = inputs['reference_fasta']
        num_threads = params["num_threads"]
        qualified_quality_phred = params["qualified_quality_phred"]
        reads_to_process = params["reads_to_process"]
        min_identity = params["min_identity"]
        min_coverage = params["min_coverage"]
        skip_read_threshold = params["skip_read_threshold"]
        min_read_length = params["min_read_length"]

        # Execute the command
        shell_proxy = FilteringShellProxyHelper.create_proxy(self.message_dispatcher)
        fastp_output_folder_path = os.path.join(shell_proxy.working_dir, 'fastp_output_folder')
        os.makedirs(fastp_output_folder_path)

        identity_coverage_output_folder_path = os.path.join(shell_proxy.working_dir, 'identity_coverage_output_folder')
        os.makedirs(identity_coverage_output_folder_path)

        final_output_folder_path = os.path.join(shell_proxy.working_dir, 'final_output_folder')
        os.makedirs(final_output_folder_path)

        summary_output_folder_path = os.path.join(shell_proxy.working_dir, 'summary_output_folder')
        os.makedirs(summary_output_folder_path)
        # call python file
        cmd = f"python3 {self.python_file_path} {input_folder.path} {fastp_output_folder_path} {identity_coverage_output_folder_path} {final_output_folder_path} {summary_output_folder_path} {int(num_threads)} {int(qualified_quality_phred)} {int(reads_to_process)} {int(min_identity)} {int(min_coverage)} {reference_fasta.path} {int(skip_read_threshold)} {int(min_read_length)}"
        shell_proxy.run(cmd, shell_mode=True)
        self.shell_proxy = shell_proxy


        folder_fastp_output = Folder(fastp_output_folder_path)
        folder_identity_coverage_output = Folder(identity_coverage_output_folder_path)
        folder_final_output = Folder(final_output_folder_path)
        folder_summary_output = Folder(summary_output_folder_path)

        summary_table_path = os.path.join(summary_output_folder_path, 'filtering_summary.csv')
        summary_table = TableImporter.call(File(summary_table_path), {'delimiter': ',', 'header': 0, 'file_format': 'csv'})

        # return the output table
        return {'final_filtered_reads': folder_final_output,
                'number_reads_filtered':summary_table
        }

    def run_after_task(self):
        if self.shell_proxy:
            self.shell_proxy.clean_working_dir()