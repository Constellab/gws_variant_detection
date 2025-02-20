# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, ConfigSpecs, File, InputSpec, Folder, ShellProxy,
                      InputSpecs, IntParam, OutputSpec, OutputSpecs, Task, TableExporter,
                      TaskInputs, TaskOutputs, task_decorator, StrParam, Table , TableImporter , ResourceSet, TaskFileDownloader)


from ..base_env.Haplotype_env import HaplotypeShellProxyHelper


@task_decorator("VariantsHaplotypesUcine", human_name="VariantsHaplotypesUcine",
                short_description="Haplotype Variant Calling")
class VariantsHaplotypesUcine(Task):
    """
   This task enables systematic matrix construction of co-occurring variants.
    """

    input_specs = InputSpecs({
        'read_variants_table':
        InputSpec(
            Table, human_name="A table containing the known variants detected per read and per barcode"),
    })

    output_specs = OutputSpecs({
        'counts_matrices_folder':
        OutputSpec(
            Folder, human_name="co-occurence matrix in term of counts",
            short_description="The frequency of reads of a target gene with multiple variants (mutation1 + mutations 2) in term of COUNTS"),

        'frequencies_matrices_folder':
        OutputSpec(
            Folder, human_name="co-occurence matrix in term of frequencies",
            short_description="The frequency of reads of a target gene with multiple variants (mutation1 + mutations 2) in term of FREQUENCY"),
        'co-occurence_table_folder':
        OutputSpec(
            Folder, human_name="co-occurence table variants",
            short_description="The frequency and count of variants across a barcode"),
        'resource_set_COUNT':
        OutputSpec(
            ResourceSet, human_name="Count Co-occurence matrix",
            short_description="Multiple Variant Detection matrices in term of COUNTS"),
        'resource_set_FREQ':
        OutputSpec(
            ResourceSet, human_name="FREQ Co-occurence matrix",
            short_description="Multiple Variant Detection matrices in term of FREQUENCIES"),
        'resource_set_cooccurence_variant':
        OutputSpec(
            ResourceSet, human_name="Co-occurence variants table",
            short_description="Multiple Variant Detection Tables in term of FREQUENCIES and COUNTS"),
    })

    config_specs: ConfigSpecs = {
    }

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_variant_haplotypes_ucine.py"
    )

    shell_proxy: ShellProxy = None


    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        # retrive the input table
        read_variants_table: Table = inputs['read_variants_table']

        # Execute the command
        shell_proxy = HaplotypeShellProxyHelper.create_proxy(self.message_dispatcher)

        counts_matrices_folder_path = os.path.join(shell_proxy.working_dir, 'counts_matrices_folder')
        os.makedirs(counts_matrices_folder_path)

        frequencies_matrices_folder_path = os.path.join(shell_proxy.working_dir, 'frequencies_matrices_folder')
        os.makedirs(frequencies_matrices_folder_path)

        co_occurence_variants_folder_path = os.path.join(shell_proxy.working_dir, 'co_occurence_variants_folder')
        os.makedirs(co_occurence_variants_folder_path)

        read_variants_file: File = TableExporter.call(read_variants_table, {
            'file_format': 'csv',
            'delimiter' : ',',
            'write_header': True,
            'write_index' : False
        })


        # call python file
        cmd = f"python3 {self.python_file_path} {read_variants_file.path} {counts_matrices_folder_path} {frequencies_matrices_folder_path} {co_occurence_variants_folder_path}"
        result = shell_proxy.run(cmd, shell_mode=True)

        if result != 0:
            raise Exception("Error during command, please check logs")
        self.shell_proxy = shell_proxy

        # display outputs :
        counts_matrices_folder = Folder(counts_matrices_folder_path)
        frequencies_matrices_folder = Folder(frequencies_matrices_folder_path)
        co_occurence_variants_folder = Folder(co_occurence_variants_folder_path)

        ## create a ressource_set
        resource_set_COUNT = ResourceSet()
        for filename in os.listdir(counts_matrices_folder_path):
            file_path = os.path.join(counts_matrices_folder_path, filename)
            if os.path.isfile(file_path) and filename.endswith('.csv'):
                # Use TableImporter.call to read CSV files and add them to the resource set
                table = TableImporter.call(File(file_path))
                resource_set_COUNT.add_resource(table, filename)

        ## create a ressource_set_freq
        resource_set_FREQ = ResourceSet()
        for filename in os.listdir(frequencies_matrices_folder_path):
            file_path = os.path.join(frequencies_matrices_folder_path, filename)
            if os.path.isfile(file_path) and filename.endswith('.csv'):
                # Use TableImporter.call to read CSV files and add them to the resource set
                table = TableImporter.call(File(file_path))
                resource_set_FREQ.add_resource(table, filename)

        resource_set_cooccurence_variant = ResourceSet()
        for filename in os.listdir(co_occurence_variants_folder_path):
            file_path = os.path.join(co_occurence_variants_folder_path, filename)
            if os.path.isfile(file_path) and filename.endswith('.csv'):
                # Use TableImporter.call to read CSV files and add them to the resource set
                table = TableImporter.call(File(file_path))
                resource_set_cooccurence_variant.add_resource(table, filename)

        # Return the outputs
        return {
            'counts_matrices_folder':counts_matrices_folder,
            'frequencies_matrices_folder':frequencies_matrices_folder,
            'co-occurence_table_folder':co_occurence_variants_folder,
            'resource_set_COUNT':resource_set_COUNT,
            'resource_set_FREQ':resource_set_FREQ,
            'resource_set_cooccurence_variant':resource_set_cooccurence_variant,
        }

    def run_after_task(self):
        if self.shell_proxy:
            self.shell_proxy.clean_working_dir()