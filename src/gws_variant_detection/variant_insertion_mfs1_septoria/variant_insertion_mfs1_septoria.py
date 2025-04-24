# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, ConfigSpecs, File, InputSpec, Folder,ShellProxy,
                      InputSpecs, IntParam, OutputSpec, OutputSpecs, Task,
                      TaskInputs, TaskOutputs, task_decorator, StrParam, Table , TableImporter, StrParam, FloatParam)


from ..base_env.Mfs_env import MfsShellProxyHelper


@task_decorator("SeptoriaInsertionsMfs1", human_name="SeptoriaInsertionsMfs1",
                short_description="Determination of insertion size of the MFS1 gene in Septoria.")

class SeptoriaInsertionsMfs1(Task):
    """
    This script classifies reads from filtered fastq.gz files based on their insertion sizes relative to a reference genome and generates summary statistics table. It calculates insertion sizes, classifies reads into defined types, and outputs the results in tsv format
    """

    input_specs = InputSpecs({
        'input_folder':
        InputSpec(
            Folder, human_name="Folder containing the fastq.gz files"),
        'reference_fasta':
        InputSpec(
            File, human_name="Combined reference sequences in FASTA format")
    })

    output_specs = OutputSpecs({
        'classification_summary': OutputSpec(Table, human_name="Classification summary", short_description="Table containing a summary of the read classifications for each sample"),
    })


    config_specs: ConfigSpecs = {
        "threads": IntParam(default_value=10, min_value=4, short_description=" threads"),
        "min_cov": IntParam(default_value=50, min_value=50, short_description="min coverage")
    }

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_variant_insertion_mfs1_septoria.py"
    )

    shell_proxy: ShellProxy = None

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        # retrive the input table
        input_folder: Folder = inputs['input_folder']
        reference_fasta: File = inputs['reference_fasta']
        threads = params["threads"]
        min_cov = params["min_cov"]

        # Execute the command
        shell_proxy = MfsShellProxyHelper.create_proxy(self.message_dispatcher)
        out_dir_path = os.path.join(shell_proxy.working_dir, 'out_dir_files')
        os.makedirs(out_dir_path)


        # call python file
        cmd = f"python3 {self.python_file_path} {input_folder.path} {reference_fasta.path} {out_dir_path} {int(threads)} {int(min_cov)}"
        shell_proxy.run(cmd, shell_mode=True)
        self.shell_proxy = shell_proxy

        out_dir_files = Folder(out_dir_path)

        # Read the known and unknown variants from the output path
        classification_summary_file_path = os.path.join(out_dir_path, 'merged_summary.tsv')

        classification_summary_table = TableImporter.call(File(classification_summary_file_path), {'delimiter': 'tab', 'header': 0, 'file_format': 'tsv', 'index_column': 0})

        # Return the outputs
        return {
            'classification_summary': classification_summary_table,

        }

    def run_after_task(self):
        if self.shell_proxy:
            self.shell_proxy.clean_working_dir()