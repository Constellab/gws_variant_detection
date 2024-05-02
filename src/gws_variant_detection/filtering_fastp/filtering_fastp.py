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
                short_description="#######################")
class NanoQualityCheck(Task):
    """
    This task uses Fastp
    """

    input_specs = InputSpecs({
        'input_folder':
        InputSpec(
            Folder, human_name="Folder containing fastq files")
    })

    output_specs = OutputSpecs({
        'output_path':
        OutputSpec(
            Folder, human_name="filtered reads",
            short_description="#################"),
        'filtered_output_path':
        OutputSpec(
            Folder, human_name="randomly filtered reads",
            short_description="#################"),
    })

    config_specs: ConfigSpecs = {
        "num_threads": IntParam(default_value=2, min_value=2, short_description="Number of threads "),
        "qualified_quality_phred": IntParam(default_value=10, min_value=4, short_description="qualified quality phred "),
        "reads_to_process": IntParam(default_value=5000, min_value=1000, short_description="number of reads to process")
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

        # Execute the command
        shell_proxy = VariantsShellProxyHelper.create_proxy(self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, 'Filtered_reads')
        os.makedirs(result_path)

        filtered_result_path = os.path.join(shell_proxy.working_dir, 'Filtered_reads_test')
        os.makedirs(filtered_result_path)
        # call python file
        cmd = f"python3 {self.python_file_path} {input_folder.path} {result_path} {filtered_result_path} {int(num_threads)} {int(qualified_quality_phred)} {int(reads_to_process)} "
        shell_proxy.run(cmd, shell_mode=True)


        folder = Folder(result_path)
        folder_filtered = Folder(filtered_result_path)

        # return the output table
        return {'output_path': folder,
                'filtered_output_path': folder_filtered}
