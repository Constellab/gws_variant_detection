# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, ConfigSpecs, File, InputSpec, Folder, ShellProxy,
                      InputSpecs, OutputSpec, OutputSpecs, Task,
                      TaskInputs, TaskOutputs, task_decorator)


from ..base_env.Mappy_env import MappyShellProxyHelper


@task_decorator("ReadsClassification", human_name="ReadsClassification",
                short_description="Per gene read Classification")
class RawReadsClassification(Task):
    """
    This task permit to classify sequencing reads based on their alignment to reference sequences, leveraging the Python bindings for Minimap2 provided by Mappy.
    The script operates by processing compressed FASTQ files located in a specified input directory, aligning each read to a reference sequence using an identity threshold to determine match quality.
    If a read aligns to the reference sequence with an identity percentage equal or higher than the set threshold, it is classified under the gene corresponding to that reference sequence.
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
            short_description="Each folder contains the different reads that map to a specific reference sequence."),
    })

    config_specs = ConfigSpecs({})

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_reads_classification.py"
    )

    shell_proxy: ShellProxy = None

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
        self.shell_proxy = shell_proxy


        folder = Folder(result_path)

        # return the output table
        return {'output_path': folder}

    def run_after_task(self):
        if self.shell_proxy:
            self.shell_proxy.clean_working_dir()