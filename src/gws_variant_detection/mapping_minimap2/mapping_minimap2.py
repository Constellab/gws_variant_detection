# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, ConfigSpecs, File, InputSpec, Folder,
                      InputSpecs, IntParam, OutputSpec, OutputSpecs, Task,
                      TaskInputs, TaskOutputs, task_decorator)


from ..base_env.Variants_env import VariantsShellProxyHelper


@task_decorator("ReadsMapping", human_name="NanoReadsMapping",
                short_description="Mapping of long reads to reference sequence")
class NanoMapping(Task):
    """
    This task is dedicated to mapping sequencing reads against a reference genome, a crucial step in the variant calling pipeline using Minimap2.
    Minimap2 is optimized for aligning long reads, which can be challenging due to their high error rate and complexity.
    """

    input_specs = InputSpecs({
        'input_folder':
        InputSpec(
            Folder, human_name="Folder containing fastq files"),
        'input_file':
        InputSpec(
            File, human_name="File containing reference sequences")
    })

    output_specs = OutputSpecs({
        'output_path':
        OutputSpec(
            Folder, human_name="Aligned reads",
            short_description="alignment files in BAM format")
    })

    config_specs: ConfigSpecs = {
        "num_threads": IntParam(default_value=2, min_value=2, short_description="Number of threads "),
    }

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_mapping_minimap2.py"
    )

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        # retrive the input table
        input_folder: Folder = inputs['input_folder']
        input_file: File = inputs['input_file']
        num_threads = params["num_threads"]

        # Execute the command
        shell_proxy = VariantsShellProxyHelper.create_proxy(self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, 'Aligned_reads')
        os.makedirs(result_path)

        # call python file
        cmd = f"python3 {self.python_file_path} {input_folder.path} {input_file.path } {int(num_threads)} {result_path}"
        shell_proxy.run(cmd, shell_mode=True)


        folder = Folder(result_path)

        # return the output table
        return {'output_path': folder}
