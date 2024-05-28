# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os

from gws_core import (ConfigParams, ConfigSpecs, File, InputSpec, Folder,
                      InputSpecs, IntParam, OutputSpec, OutputSpecs, Task,
                      TaskInputs, TaskOutputs, task_decorator, StrParam, Table , TableImporter , ResourceSet)


from ..base_env.Variants_env import VariantsShellProxyHelper


@task_decorator("VariantsHaplotypes", human_name="VariantsHaplotypes",
                short_description="Haplotype Variant Calling")
class VariantsHaplotypes(Task):
    """
   This task enables systematic matrix construction from phased variant data using Whatshap, facilitating the uncovering of intricate genetic variation patterns and thereby advancing the understanding of disease and trait genetics.
    """

    input_specs = InputSpecs({
        'input_vcf_folder':
        InputSpec(
            Folder, human_name="Folder containing the different normalized variant calling files"),
        'input_bam_folder':
        InputSpec(
            Folder, human_name="Folder containing the different alignment files"),
        'reference_genome':
        InputSpec(
            File, human_name="File containing gene reference sequence")
    })

    output_specs = OutputSpecs({
        'output_folder':
        OutputSpec(
            Folder, human_name="phased variants",
            short_description="A folder containing the different phased vcf files"),
        'output_matrix_folder':
        OutputSpec(
            Folder, human_name="Tables for detecting multiple variants",
            short_description="The frequency of reads of a target gene with multiple variants (mutation1 + mutations 2)"),
        'resource_set':
        OutputSpec(
            ResourceSet, human_name="Metadata associated to each cell",
            short_description="This table stores metadata associated with each cell such as mitochondrial content , number of counts etc"),

    })

    config_specs: ConfigSpecs = {
        "gene_name":
        StrParam(allowed_values=["Cytb","OSBP","Tubulin","ATPase"], default_value="Cytb",
                 short_description="Targeted gene")
    }

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_variant_haplotypes.py"
    )

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        """ Run the task """

        # retrive the input table
        input_vcf_folder: Folder = inputs['input_vcf_folder']
        input_bam_folder: Folder = inputs['input_bam_folder']
        reference_genome: File = inputs['reference_genome']
        gene_name = params["gene_name"]

        # Execute the command
        shell_proxy = VariantsShellProxyHelper.create_proxy(self.message_dispatcher)
        known_variants_vcf_files_path = os.path.join(shell_proxy.working_dir, 'known_variants_vcf_files')
        os.makedirs(known_variants_vcf_files_path)


        phased_known_variants_vcf_files_path = os.path.join(shell_proxy.working_dir, 'phased_known_variants_vcf_files')
        os.makedirs(phased_known_variants_vcf_files_path)

        output_matrix_folder_path = os.path.join(shell_proxy.working_dir, 'output_matrix_folder')
        os.makedirs(output_matrix_folder_path)

        # call python file
        cmd = f"python3 {self.python_file_path} {input_vcf_folder.path} {input_bam_folder.path} {reference_genome.path } {str(gene_name)} {known_variants_vcf_files_path} {phased_known_variants_vcf_files_path} {output_matrix_folder_path}"
        shell_proxy.run(cmd, shell_mode=True)

        # display outputs :
        # 1_normalized vcf
        phased_known_variants_vcf_files = Folder(phased_known_variants_vcf_files_path)
        output_matrix_folder = Folder(output_matrix_folder_path)

        ## create a ressource_set
        resource_set = ResourceSet()
        for filename in os.listdir(output_matrix_folder_path):
            file_path = os.path.join(output_matrix_folder_path, filename)
            if os.path.isfile(file_path) and filename.endswith('.csv'):
                # Use TableImporter.call to read CSV files and add them to the resource set
                table = TableImporter.call(File(file_path))
                resource_set.add_resource(table, filename)

        # Return the outputs
        return {
            'output_folder': phased_known_variants_vcf_files,
            'output_matrix_folder': output_matrix_folder,
            'resource_set':resource_set
        }