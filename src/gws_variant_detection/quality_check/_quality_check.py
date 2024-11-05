
import sys
import pandas as pd
import os
import subprocess

class NanoPlotQualityChecker:
    def __init__(self, input_folder, output_folder, num_threads):
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.num_threads = num_threads

    def check_quality(self):
        # Create the output folder if it doesn't exist
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        # Get a list of all FASTQ files in the input folder
        fastq_files = [f for f in os.listdir(self.input_folder) if f.endswith('.fastq') or f.endswith('.fastq.gz')]

        # Run NanoPlot for each FASTQ file
        for fastq_file in fastq_files:
            input_path = os.path.join(self.input_folder, fastq_file)
            output_path = os.path.join(self.output_folder, fastq_file.split('.')[0])  # Output path construction
            nanoplot_cmd = f'NanoPlot --fastq {input_path} -o {output_path} --threads {self.num_threads} --tsv_stats --raw'
            subprocess.run(nanoplot_cmd, shell=True)

class NanoPlotSummary:
    def __init__(self, quality_output_folder):
        self.quality_output_folder = quality_output_folder
        self.desired_columns = ['number_of_reads', 'number_of_bases', 'median_read_length', 'mean_read_length',
                                'Reads >Q5:', 'Reads >Q7:', 'Reads >Q10:', 'Reads >Q12:', 'Reads >Q15:']

    def clean_column_names(self, df):
        # Remove the colon ":" from specific column names
        df.columns = df.columns.str.replace(':', '', regex=True)

    def merge_nano_stats(self):
        # Initialize an empty DataFrame to store merged data
        merged_data = pd.DataFrame()

        # Recursively search for NanoStats.txt files within the quality output folder
        for root, dirs, files in os.walk(self.quality_output_folder):
            for file in files:
                if file == 'NanoStats.txt':
                    # Get the full path to the NanoStats.txt file
                    nano_stats_file = os.path.join(root, file)

                    # Read NanoStats.txt file and merge the data into the DataFrame
                    with open(nano_stats_file, 'r') as f:
                        lines = f.readlines()
                        data = {}
                        for line in lines[1:]:
                            parts = line.strip().split('\t')
                            metric = parts[0]
                            value = parts[1].split(')')[0].strip()  # Remove everything after the closing parenthesis
                            if metric in self.desired_columns:
                                if metric in ['Reads >Q5:', 'Reads >Q7:', 'Reads >Q10:', 'Reads >Q12:', 'Reads >Q15:']:
                                    value += ')'  # Preserve the closing parenthesis for these metrics
                                data[metric] = value
                        # Use the folder name as the column name
                        folder_name = os.path.basename(os.path.dirname(nano_stats_file))
                        # Convert data to DataFrame with a single column and assign it to the column_name
                        if data:
                            merged_data[folder_name] = pd.DataFrame(data.values(), index=data.keys())

        # Transpose the DataFrame to have barcodes as rows and metrics as columns
        merged_data = merged_data.transpose()

        # Extract the barcode number from the index, sort, and then reindex the DataFrame
        merged_data['Barcode_Number'] = merged_data.index.to_series().str.extract(r'barcode(\d+)-')[0].astype(int)
        merged_data.sort_values(by='Barcode_Number', inplace=True)
        merged_data.drop(columns=['Barcode_Number'], inplace=True)

        # Clean column names
        self.clean_column_names(merged_data)
        # Write the merged DataFrame to a single NanoStats.txt file
        output_file = 'NanoStats_Summary.txt'
        
        merged_data.to_csv(output_file, sep='\t')



# Example usage
input_folder = sys.argv[1]
num_threads = int(sys.argv[2])
output_folder = sys.argv[3]

nanoplot_checker = NanoPlotQualityChecker(input_folder, output_folder, num_threads)
nanoplot_checker.check_quality()

# Merge NanoStats.txt files into a single table
summary_generator = NanoPlotSummary(output_folder)
summary_generator.merge_nano_stats()

