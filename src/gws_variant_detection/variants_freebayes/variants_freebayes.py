# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, ConfigSpecs, File, InputSpec, Folder,
                      InputSpecs, IntParam, OutputSpec, OutputSpecs, Task,
                      TaskInputs, TaskOutputs, task_decorator , TableImporter , Table , FloatParam)


from ..base_env.Variants_env import VariantsShellProxyHelper


@task_decorator("FreebayesVariantsDetection", human_name="NanoFreebayesVariantsDetection",
                short_description="#######################")
class FreebayesVariantsDetection(Task):
    """
    This task uses freebayes
    ploidy refers to the number of copies of each chromosome present in the genome of the organism being analyzed.
    For diploid organisms like humans, the ploidy is 2, indicating that there are two copies of each chromosome
    (one inherited from each parent). In polyploid organisms or organisms with different levels of ploidy
    (e.g., triploid, tetraploid), the ploidy value would be different accordingly.
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
            File, human_name="Variants calling",
            short_description="#################"),
    })

    config_specs: ConfigSpecs = {
        "ploidy": IntParam(default_value=2, min_value=1, short_description="Set the  the number of copies of each chromosome present in the genome of the organism being analyzed"),
        "detection_threshold": FloatParam(default_value=0.05, min_value=0, short_description="Set the detection threshold value"),

    }

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_variants_freebayes.py"
    )

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        # retrive the input table
        input_folder: Folder = inputs['input_folder']
        input_file: File = inputs['input_file']
        ploidy = params["ploidy"]
        detection_threshold = params["detection_threshold"]

        # Execute the command
        shell_proxy = VariantsShellProxyHelper.create_proxy(self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, 'Raw_Variants_Calling')
        os.makedirs(result_path)

        # call python file
        cmd = f"python3 {self.python_file_path} {input_folder.path} {input_file.path } {result_path} {int(ploidy)} {float(detection_threshold)}"
        shell_proxy.run(cmd, shell_mode=True)


        vcf_file = os.path.join(result_path, "variants_raw.vcf")



        # return the output table
        return {
                'output_path': File(vcf_file)
        }
