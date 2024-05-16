# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, ConfigSpecs, File, InputSpec, Folder,
                      InputSpecs, IntParam, OutputSpec, OutputSpecs, Task,
                      TaskInputs, TaskOutputs, task_decorator , TableImporter , Table , FloatParam)


from ..base_env.Variants_env import VariantsShellProxyHelper


@task_decorator("MinipileupVariantsDetection", human_name="NanoMinipileupVariantsDetection",
                short_description="#######################")
class MinipileupVariantsDetection(Task):
    """
    This task uses minipileup
    """

    input_specs = InputSpecs({
        'input_folder':
        InputSpec(
            Folder, human_name="Folder containing alignments files"),
        'input_file':
        InputSpec(
            File, human_name="File containing reference sequences")
    })

    output_specs = OutputSpecs({
        'output_path':
        OutputSpec(
            Folder, human_name="Variants calling",
            short_description="#################"),
    })

    config_specs: ConfigSpecs = {
    }

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_variants_minipileup.py"
    )

    minipileup_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "minipileup"
    )

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        # retrive the input table
        input_folder: Folder = inputs['input_folder']
        input_file: File = inputs['input_file']

        # Execute the command
        shell_proxy = VariantsShellProxyHelper.create_proxy(self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, 'Raw_Variants_Calling')
        os.makedirs(result_path)

        # call python file
        cmd = f"python3 {self.python_file_path} {self.minipileup_file_path} {input_folder.path} {input_file.path} {result_path}"
        shell_proxy.run(cmd, shell_mode=True)


        vcf_file = os.path.join(result_path)



        # return the output table
        return {
                'output_path': Folder(vcf_file)
        }
