# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, ConfigSpecs, File, InputSpec, Folder, ShellProxy,
                      InputSpecs, IntParam, OutputSpec, OutputSpecs, Task,
                      TaskInputs, TaskOutputs, task_decorator , Table , TableImporter)


from ..base_env.Variants_env import VariantsShellProxyHelper


@task_decorator("QualityCheck", human_name="NanoQualitycheck",
                short_description="Statistical summary of reads quality")
class NanoQualityCheck(Task):
    """
    This task uses Nanoplot in order to check the quality of nanopore reads. It takes a folder containing fastq reads and it generates
    a table with different statistical values :
     <b style="margin-left: Number of Reads: </b> This refers to the total count of individual sequencing reads obtained from a sequencing experiment.
    Each read represents a fragment of DNA or RNA that has been sequenced.
     <b style="margin-left: Number of Bases: </b> This indicates the total count of nucleotide bases (A, T, C, G, and sometimes N) across all sequencing reads.
    It represents the cumulative length of all the sequences obtained in a sequencing experiment.
     <b style="margin-left: </b> Median Read Length: The median read length is the middle value of all read lengths when they are arranged in ascending order.
    It provides a measure of the typical or central read length in a dataset.
     <b style="margin-left: Mean Read Length: </b> The mean read length is the average length of all sequencing reads in a dataset.
    It is calculated by dividing the total number of bases by the number of reads.
     <b style="margin-left: Reads >Q10 , Reads >Q15 , Reads >Q20: </b> This refers to the number of reads with a quality score respectively greater than 5 , 7 , 10 , 12 and 15.
    Quality scores represent the confidence level of each base call in a sequencing read.
    """

    input_specs = InputSpecs({
        'input_folder':
        InputSpec(
            Folder, human_name="Folder containing fastq files")
    })

    output_specs = OutputSpecs({
        'output_path':
        OutputSpec(
            Table, human_name="quality check statistical table",
            short_description="This table contains the different statistical values")
    })

    config_specs: ConfigSpecs = {
        "num_threads": IntParam(default_value=15, min_value=2, short_description="Number of threads ")
    }

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_quality_check.py"
    )

    shell_proxy: ShellProxy = None

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        # retrive the input table
        input_folder: Folder = inputs['input_folder']
        num_threads = params["num_threads"]

        # Execute the command
        shell_proxy = VariantsShellProxyHelper.create_proxy(self.message_dispatcher)
        result_path = os.path.join(shell_proxy.working_dir, 'quality_plots')
        os.makedirs(result_path)

        # call python file
        cmd = f"python3 {self.python_file_path} {input_folder.path} {int(num_threads)} {result_path} "
        shell_proxy.run(cmd, shell_mode=True)
        self.shell_proxy = shell_proxy

        statistical_file_name = os.path.join(
            shell_proxy.working_dir, "NanoStats_Summary.txt")

        statistical_table = TableImporter.call(
            File(statistical_file_name),
            {'delimiter': 'tab', 'header': 0, 'index_column': 0})


        folder = Folder(result_path)

        # return the output table
        return {'output_path': statistical_table}

    def run_after_task(self):
        if self.shell_proxy:
            self.shell_proxy.clean_working_dir()