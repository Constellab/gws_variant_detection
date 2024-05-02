# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, ConfigSpecs, File, InputSpec, Folder,
                      InputSpecs, OutputSpec, OutputSpecs, Task,
                      TaskInputs, TaskOutputs, task_decorator)


from ..base_env.Mappy_env import MappyShellProxyHelper


@task_decorator("ReadsClassification", human_name="ReadsClassification",
                short_description="Classify each read")
class RawReadsClassification(Task):
    """
    This task uses mappy the Python binding for minimap2, which allows you to use minimap2's alignment functionality directly from Python.
    This is what the script uses to align reads to reference sequences.
    """

    input_specs = InputSpecs({
        'input_folder':
        InputSpec(
            Folder, human_name="Folder containing different barcodes folder"),
        'input_file':
        InputSpec(
            File, human_name="Merged reference sequences")
    })

    output_specs = OutputSpecs({
        'output_path':
        OutputSpec(
            Folder, human_name="per gene folder",
            short_description="#################"),
    })

    config_specs: ConfigSpecs = {
    }

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_reads_classification.py"
    )

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        # retrive the input table
        input_folder: Folder = inputs['input_folder']
        input_file: File = inputs['input_file']

        # Execute the command
        shell_proxy = MappyShellProxyHelper.create_proxy(self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, 'Raw_Reads_Classification')
        os.makedirs(result_path)

        # call python file
        cmd = f"python3 {self.python_file_path} {input_folder.path} {input_file.path } {result_path}"
        shell_proxy.run(cmd, shell_mode=True)


        folder = Folder(result_path)

        # return the output table
        return {'output_path': folder}