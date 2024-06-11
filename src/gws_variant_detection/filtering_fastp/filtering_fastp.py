# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, ConfigSpecs, File, InputSpec, Folder,
                      InputSpecs, IntParam, OutputSpec, OutputSpecs, Task,
                      TaskInputs, TaskOutputs, task_decorator)


from ..base_env.Variants_env import VariantsShellProxyHelper


@task_decorator("ReadsQualityFiltering", human_name="NanoReadsFiltering",
                short_description="Elimination of low-quality reads and fixing the number of reads per sample.")
class ReadsQualityFiltering(Task):
    """
    This task automates the process of filtering and randomly sampling reads from sequencing data using fastp and seqtk respectively.
    The number of threads, a quality threshold, and the number of reads to process must be specified.
    """

    input_specs = InputSpecs({
        'input_folder':
        InputSpec(
            Folder, human_name="Folder containing fastq files")
    })

    output_specs = OutputSpecs({
        'output_path':
        OutputSpec(
            Folder, human_name="Filtered reads",
            short_description="This folder contains raw filtered reads having quality > threshold"),
        'filtered_output_path':
        OutputSpec(
            Folder, human_name="Subset of filtered reads",
            short_description="This folder contains a subset of randomly filtered reads"),
    })

    config_specs: ConfigSpecs = {
        "num_threads": IntParam(default_value=15, min_value=2, short_description="Number of threads "),
        "qualified_quality_phred": IntParam(default_value=10, min_value=4, short_description="qualified quality phred "),
        "reads_to_process": IntParam(default_value=5000, min_value=1000, short_description="number of reads to process"),
        "seqtk_sample_seed": IntParam( default_value=100, min_value=1, human_name="sample_seed", visibility=IntParam.PROTECTED_VISIBILITY,
            short_description="Fix sample seed")
    }

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_filtering_fastp.py"
    )

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        # retrive the input table
        input_folder: Folder = inputs['input_folder']
        num_threads = params["num_threads"]
        qualified_quality_phred = params["qualified_quality_phred"]
        reads_to_process = params["reads_to_process"]
        seqtk_sample_seed = params["seqtk_sample_seed"]

        # Execute the command
        shell_proxy = VariantsShellProxyHelper.create_proxy(self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, 'Filtered_reads')
        os.makedirs(result_path)

        filtered_result_path = os.path.join(shell_proxy.working_dir, 'fixed_number_of_filtered_reads')
        os.makedirs(filtered_result_path)
        # call python file
        cmd = f"python3 {self.python_file_path} {input_folder.path} {result_path} {filtered_result_path} {int(num_threads)} {int(qualified_quality_phred)} {int(reads_to_process)} {int(seqtk_sample_seed)}"
        shell_proxy.run(cmd, shell_mode=True)


        folder = Folder(result_path)
        folder_filtered = Folder(filtered_result_path)

        # return the output table
        return {'output_path': folder,
                'filtered_output_path': folder_filtered}
