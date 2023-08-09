# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

import requests
import sys
import pandas as pd
from math import ceil
import os
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s:.%(funcName)s: %(message)s')
logger = logging.getLogger(__name__)

#import argparse

#PARAMS
project: str = "TCGA-LIHC"
strategy: str = "WXS"
format: str = "BAM"
slice_size: int = 200
#END PARAMS


def retreiveProjectData(project: str, strategy: str, format: str, output_name: str, manifest: bool = False) -> None:
    #Format: BAM, MAF, expression
    #project: Any TCGA project ID e.g. TCGA-LIHC
    #strategy: WXS, WGS, RNA-Seq, miRNA-Seq, etc.

    fields = [
        "file_name",
        "total_reads",
        "cases.submitter_id",
        "average_read_length",
        "mean_coverage",
        "pairs_on_different_chr",
        "proportion_targets_no_coverage",
        "proportion_reads_mapped",

        "cases.case_id",
        "cases.disease_type",

        "cases.project.project_id",
        "cases.project.project_name",

        "cases.demographic.vital_status",
        "cases.demographic.days_to_birth",
        "cases.demographic.days_to_death",

        "cases.samples.sample_type",
        "cases.samples.days_to_collection",
        "cases.samples.days_to_sample_procurment",
        "cases.samples.portions.analytes.aliquots.concentration",
        "cases.samples.portions.analytes.concentration",
        "cases.samples.sample_type_id",
        "cases.samples.submitter_id",
        "cases.samples.preservation_method",
        "cases.samples.tissue_type",
        "cases.samples.tumor_descriptor",
        "cases.samples.specimen_type",

        "cases.diagnoses.days_to_last_followup",
        "cases.diagnoses.tumor_grade",
        "cases.diagnoses.tumor_stage",
        "cases.diagnoses.age_at_diagnosis",
        "cases.diagnoses.ajcc_pathologic_stage",
        "cases.diagnoses.ajcc_pathologic_n",
        "cases.diagnoses.ajcc_pathologic_m",
        "cases.diagnoses.ajcc_pathologic_t",
        "cases.diagnoses.figo_stage",
        "cases.diagnoses.prior_malignancy",
        "cases.diagnoses.prior_treatment",
        ]

    fields = ",".join(fields)

    files_endpt = "https://api.gdc.cancer.gov/files" #files

    # This set of filters is nested under an 'and' operator.
    if format=="BAM":

        filters = {
            "op": "and",
            "content":[
                {
                "op": "in",
                "content":{
                    "field": "cases.project.project_id",
                    "value": [project]
                    }
                },
                {
                "op": "in",
                "content":{
                    "field": "files.experimental_strategy",
                    "value": [strategy]
                    }
                },
                {
                "op": "in",
                "content":{
                    "field": "files.data_format",
                    "value": [format]
                    }
                }
            ]
        }

    elif format=="MAF":

        filters = {
            "op": "and",
            "content":[
                {
                "op": "in",
                "content":{
                    "field": "cases.project.project_id",
                    "value": [project]
                    }
                },
                {
                "op": "in",
                "content":{
                    "field": "files.experimental_strategy",
                    "value": [strategy]
                    }
                },
                {
                "op": "in",
                "content":{
                    "field": "files.data_format",
                    "value": [format]
                    }
                },
                {
                "op": "in",
                "content":{
                    "field": "files.data_type",
                    "value": "Masked Somatic Mutation"
                }
            }
        ]
    }

    elif format=="expression":
        
        filters = {
            "op": "and",
            "content":[
                {
                "op": "in",
                "content":{
                    "field": "cases.project.project_id",
                    "value": [project]
                    }
                },
                {
                "op": "in",
                "content":{
                    "field": "files.data_type",
                    "value": "Gene Expression Quantification"
                    }
                 }
            ]
        }

    else:
        logging.error("ERROR: format not recognized")
        return

    # A POST is used, so the filter parameters can be passed directly as a Dict object.
    params = {
        "filters": filters,
        "fields": fields,
        "format": "TSV",
        "size": "10000"
        }

    if manifest:
        params = {
        "filters": filters,
        "fields": fields,
        "format": "TSV",
        "size": "10000",
        "return_type":"manifest"
        }

    # The parameters are passed to 'json' rather than 'params' in this case
    response = requests.post(files_endpt, headers = {"Content-Type": "application/json"}, json = params)
    original_stdout = sys.stdout

    with open(output_name, 'w') as f:
        sys.stdout = f # Change the standard output to the file we created.
        logging.info(response.content.decode("utf-8"))
        sys.stdout = original_stdout # Reset the standard output to its original value

def createManifestSlices(samples_file_metadata: str, samples_file_manifest: str, slice_size: int, output_name_root: str) -> None:
    """
    Slices a TCGA manifest file into smaller manifest files containing a specified number of file IDs.

    Args:
        samples_file_metadata (str): The path to the TCGA metadata file.
        samples_file_manifest (str): The path to the TCGA manifest file.
        slice_size (int): The number of file IDs to include in each sliced manifest file.
        output_name_root (str): The root name for the output sliced manifest files.

    Returns:
        None: This function does not return anything, it simply creates the sliced manifest files.

    Raises:
        None: This function does not raise any exceptions.

    Example:
        createManifestSlices("metadata.txt", "manifest.txt", 100, "sliced_manifest")
    """
    #open the file with IDs load to pandas df
    df_metadata = pd.read_csv(samples_file_metadata, sep="\t")
    df_metadata = df_metadata.sort_values(by=["cases.0.case_id"])
    id_list = df_metadata["id"].tolist()

    id_list_length = len(id_list)
    slices = ceil(id_list_length/slice_size)

    logging.info(f"Found {id_list_length} file IDs. Creating {slices} manifest files with {slice_size} IDs each")

    df_manifest = pd.read_csv(samples_file_manifest, sep="\t")

    slices_dir="manifest_slices"

    if not os.path.exists(slices_dir):
        os.makedirs(slices_dir)

    for i in range(slices):
        #slice the table
        file_uuid_list = id_list[(i*slice_size):((i+1)*slice_size)]

        df_manifest_sliced = df_manifest[df_manifest["id"].isin(file_uuid_list)]

        df_manifest_sliced.to_csv(f"{slices_dir}/{output_name_root}-slice-{i}.txt", index=False, sep="\t")

    logging.info("Finished!")

def main(project: str, strategy: str, format: str, slice_size: int):
    """
    Downloads TCGA data for a specified project, strategy, and data format, and slices the manifest file into smaller files.

    Args:
        project (str): The TCGA project ID to download data for.
        strategy (str): The sequencing strategy to download data for (e.g. "RNASeq", "WGS", "WXS").
        format (str): The data format to download (e.g. "BAM", "MAF", "expression").
        slice_size (int): The number of file IDs to include in each sliced manifest file.

    Returns:
        None: This function does not return anything, it simply downloads and slices the TCGA data.

    Raises:
        None: This function does not raise any exceptions.

    Example:
        main("TCGA-BRCA", "RNASeq", "expression", 2000)
        or
        main("TCGA-BRCA", "WXS", "BAM", 200)
    """
    # retreive project metadata
    output_name_meta = f'{project}-{strategy}-{format}-metadata.tsv'
    retreiveProjectData(project, strategy, format, output_name_meta, manifest=False)

    # retreive project manifest
    output_name_manifest = f'{project}-{strategy}-{format}-manifest.tsv'
    retreiveProjectData(project, strategy, format, output_name_manifest, manifest=True)

    # slice the project manifest by ordered metadata case-ids into chunks of size slice_size
    samples_file_manifest = f"{output_name_manifest}"
    samples_file_metadata = f"{output_name_meta}"

    if format != "expression": 
        output_name_root = f"{project}-{strategy}-{format}-manifest"
        createManifestSlices(samples_file_metadata, samples_file_manifest, slice_size, output_name_root)

if __name__ == "__main__":

    '''parser = argparse.ArgumentParser(description='Retreive GDC project data')
    parser.add_argument('-p', '--project', type=str, required=True, help='GDC project ID')
    parser.add_argument('-s', '--strategy', type=str, required=True, help='GDC experimental strategy, e.g. WGS, WXS, RNA-Seq, miRNA-Seq, Bisulfite-Seq, etc.')
    parser.add_argument('-f', '--format', type=str, required=True, help='GDC data format, e.g. BAM, MAF, expression (for gene expression quantification)')
    parser.add_argument('-n', '--slice_size', type=int, required=True, help='Size of manifest slices')

    args = parser.parse_args()

    main(args.project, args.strategy, args.format, args.slice_size)'''

    main(project, strategy, format, slice_size)