# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, ConfigSpecs, File, InputSpec, Folder,
                      InputSpecs, IntParam, OutputSpec, OutputSpecs, Task,
                      TaskInputs, TaskOutputs, task_decorator , Table , TableImporter)


from ..base_env.Variants_env import VariantsShellProxyHelper


@task_decorator("PorechopTrimming", human_name="NanoAdaptersElimination",
                short_description="Adaptors elimination" , hide=True)
class PorechopTrimming(Task):
    """
    This task uses Porechop to trim reads and eliminate adapters
    """

    input_specs = InputSpecs({
        'input_folder':
        InputSpec(
            Folder, human_name="Folder containing fastq files")
    })

    output_specs = OutputSpecs({
        'output_path':
        OutputSpec(
            Folder, human_name="trimmed reads",
            short_description="#################")
    })

    config_specs = ConfigSpecs({
        "num_threads": IntParam(default_value=2, min_value=2, short_description="Number of threads ")
    })

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_trimming_porechop.py"
    )

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        # retrive the input table
        input_folder: Folder = inputs['input_folder']
        num_threads = params["num_threads"]

        # Execute the command
        shell_proxy = VariantsShellProxyHelper.create_proxy(self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, 'trimmed_reads')
        os.makedirs(result_path)

        # call python file
        cmd = f"python3 {self.python_file_path} {input_folder.path} {int(num_threads)} {result_path} "
        shell_proxy.run(cmd, shell_mode=True)


        folder = Folder(result_path)

        # return the output table
        return {'output_path': folder}
