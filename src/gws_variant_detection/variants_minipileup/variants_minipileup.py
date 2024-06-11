# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, ConfigSpecs, File, InputSpec, Folder,
                      InputSpecs, IntParam, OutputSpec, OutputSpecs, Task,
                      TaskInputs, TaskOutputs, task_decorator, TableImporter, Table, FloatParam, TaskFileDownloader)


from ..base_env.Variants_env import VariantsShellProxyHelper
import subprocess


@task_decorator("MinipileupVariantsDetection", human_name="NanoMinipileupVariantsDetection",
                short_description="Variant calling using minipileup")
class MinipileupVariantsDetection(Task):
    """
    This task uses Minipileup a simple pileup-based variant caller.
    It iteratively processes each position in the genome where reads are aligned (pileup). A 'Pileup' is when a multiple sequencing reads overlap at the same genomic position
    Minipileup counts the number of each type of base (allele) present at every genomic position. This count is crucial to determine if a variation is present and to assess the allele frequency.
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
            short_description="per base variant calling file"),
    })

    config_specs: ConfigSpecs = {
    }

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_variants_minipileup.py"
    )

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        # retrive the input table
        input_folder: Folder = inputs['input_folder']
        input_file: File = inputs['input_file']

        # Setup the minipileup executable
        downloader = TaskFileDownloader(
            brick_name=MinipileupVariantsDetection.get_brick_name(),
            message_dispatcher=self.message_dispatcher)
        minipileup_folder = downloader.download_file_if_missing(
            url='https://storage.sbg.cloud.ovh.net/v1/AUTH_a0286631d7b24afba3f3cdebed2992aa/public/gws_variant_detection/minipileup-1.0.zip',
            filename='minipileup', decompress_file=True)

        minipileup_folder_v1 = os.path.join(os.path.join(minipileup_folder, 'minipileup-1.0'))

        #  path including the directory if exists
        minipileup_file_path = os.path.join(minipileup_folder_v1, 'minipileup')

        if not os.path.isfile(minipileup_file_path):
            # If the executable does not exist, navigate to the directory and run the make command
            shell_proxy = VariantsShellProxyHelper.create_proxy(self.message_dispatcher, minipileup_folder_v1)

            # Execute the make command
            make_command = "make"
            result = shell_proxy.run(make_command, shell_mode=True)
            if result:
                raise Exception("Error during build of minipileup, please check logs")

        # Execute the command
        shell_proxy = VariantsShellProxyHelper.create_proxy(self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, 'Raw_Variants_Calling')
        os.makedirs(result_path)

        # call python file
        cmd = f"python3 {self.python_file_path} {minipileup_file_path} {input_folder.path} {input_file.path} {result_path}"
        shell_proxy.run(cmd, shell_mode=True)

        vcf_file = os.path.join(result_path)

        # return the output table
        return {
            'output_path': Folder(vcf_file)
        }
