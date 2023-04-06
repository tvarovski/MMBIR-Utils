#!/Users/twarowski/mambaforge/bin/python

# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

import pandas as pd
import sys
import cancer_config as cfg
import logging


#set logging level
logging.basicConfig(level=logging.INFO, format='%(levelname)s:.%(funcName)s: %(message)s')
#example logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(name)s.%(funcName)s +%(lineno)s: %(levelname)-8s [%(process)d] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

def help():
    '''
This is a script designed to be used from a command line which takes
a command line argument and runs a corresponding functionality from the
MMBIR-UTILS package

Usage example(from the command line):
    
    python utils.py <KEYWORD>

Currently available KEYWORDs:
    - groupCases
        # groups results of the MMBSearch snakemake pipeline into Case_ID dirs

    - parseOutputs
        # creates a consolidated_results file

    - masked_snv_mv
        # moves snvs data from /nfsscratch into Case_ID dirs

    - createFullCancerTable
        # creates annotated master_mmbir_tables for raw and filtered
        # results from all available data
        # creates a manifest file with missing bam_files

    - findCosmicGenes
        #filteres out and prints out genes in the Cosmic_DB that
        #are present in the specified DF path (sys.argv[2])

    - performDiffExprAnalysis
        #!!!add parallelization!!!
        # performs differential expression analysis

    - getMissingBams
        # creates a manifest file with missing bam_files if any are missing

    - expressionParser
        # parses the expression data from the TCGA database

    - performSNVanalysis
        # performs SNV analysis

    - heatMapper
        # creates a heatmap from the specified DF path (sys.argv[2])
    
    - help
        # prints this help message

'''
    print(help.__doc__)
    
def groupCasesInit():
    from tools import groupCases

    cancer = cfg.settings["TCGA-PROJECT"]
    username = cfg.settings["username"]
    metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
    metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"
    df_metadata = pd.read_csv(metadata_location, sep="\t")

    groupCases(df_metadata)

def parseOutputsInit():
    from tools import parseOutputs

    cancer = cfg.settings["TCGA-PROJECT"]
    username = cfg.settings["username"]
    metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
    metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"
    df_metadata = pd.read_csv(metadata_location, sep="\t")

    parseOutputs(df_metadata, consolidated_results_name="consolidated_results.tsv")

def createFullCancerTableInit():
    from tools import createFullCancerTable

    cancer = cfg.settings["TCGA-PROJECT"]
    username = cfg.settings["username"]
    metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
    metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"
    df_metadata = pd.read_csv(metadata_location, sep="\t")
    manifest_file = f"TCGA-{cancer}-WXS-BAM-manifest.tsv"

    params={
        "df_metadata": df_metadata,
        "manifest_location": f"/Users/{username}/MMBIR_Databases/TCGA/{manifest_file}",
        "missing_manifest_output_name": f"missing_bams_{cancer}_manifest.tsv",
        "output_raw_name": "raw_mmbir_results_master_df.csv",
        "output_filtered_name": "filtered_mmbir_results_master_df.csv"
        }
        
    createFullCancerTable(params)

def masked_snv_mvInit():
    from tools import masked_snv_mv

    cancer = cfg.settings["TCGA-PROJECT"]
    username = cfg.settings["username"]
    metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
    metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"
    df_metadata = pd.read_csv(metadata_location, sep="\t")
    metadata_file_maf = f"TCGA-{cancer}-WXS-MAF-metadata.tsv"
    metadata_location_maf = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file_maf}"
    df_metadata = pd.read_csv(metadata_location_maf, sep="\t")

    #need to add a sys.argv[2] for the snvs_loc
    snvs_loc = sys.argv[2]

    masked_snv_mv(df_metadata, snvs_loc=snvs_loc)

def findCosmicGenesInit():
    from tools import findCosmicGenes

    census_dir = cfg.settings['cosmicdb_dir']
    mmb_df_input_path = sys.argv[2]

    filter_dict={
        "ref_complexity_filter": cfg.settings["ref_complexity_filter"],
        "bir_complexity_filter": cfg.settings["bir_complexity_filter"],
        "ref_homology_check_filter": cfg.settings["ref_homology_check_filter"],
        "bir_homology_check_filter": cfg.settings["bir_homology_check_filter"],
        "exones_only": cfg.settings["exones_only"]
        }

    findCosmicGenes(mmb_df_input_path, filter_dict, census_dir)

def performDiffExprAnalysisInit():
    from tools import performDiffExprAnalysis

    cancer = cfg.settings["TCGA-PROJECT"]
    username = cfg.settings["username"]
    metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
    metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"
    df_metadata = pd.read_csv(metadata_location, sep="\t")

    params = {
        "cancer": cancer,
        "df_metadata": df_metadata,
        "MMBIR_THRESHOLD_LOW": cfg.settings["MMBIR_THRESHOLD_LOW"],
        "MMBIR_THRESHOLD_HIGH":cfg.settings["MMBIR_THRESHOLD_HIGH"],
        "consolidated_results_path": cfg.settings["consolidated_results_path"],
        "expression_df_path": f"expression_data_{cancer}.pickle",
        "min_concentration": cfg.settings["min_concentration"],
        "mmbir_fraction_high": cfg.settings["mmbir_fraction_high"],
        "mmbir_fraction_low": cfg.settings["mmbir_fraction_low"],
        "outputs_path": cfg.settings["outputs_path"],
        "expression_df_metadata_path": f"/Users/{username}/MMBIR_Databases/TCGA/TCGA-{cancer}-WXS-expression-metadata.tsv",
        "investigated_tissue": cfg.settings["investigated_tissue"],
            }

    performDiffExprAnalysis(params)

def getMissingBamsInit():
    from tools import getMissingBams

    cancer = cfg.settings["TCGA-PROJECT"]
    username = cfg.settings["username"]
    metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
    metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"
    df_metadata = pd.read_csv(metadata_location, sep="\t")
    manifest_file = f"TCGA-{cancer}-WXS-BAM-manifest.tsv"

    params={
        "df_metadata": df_metadata,
        "manifest_location": f"/Users/{username}/MMBIR_Databases/TCGA/{manifest_file}",
        "missing_manifest_output_name": f"missing_bams_{cancer}_manifest.tsv"
        }
        
    getMissingBams(params)

def expressionParserInit():
        
    from tools import expressionParser

    cancer = cfg.settings["TCGA-PROJECT"]
    username = cfg.settings["username"]

    expression_data_path_root = cfg.settings["expression_data_path_root"]
    expression_data_path = f"{expression_data_path_root}/{username}/TCGA-{cancer}/expression"

    expression_metadata_file = f"TCGA-{cancer}-WXS-expression-metadata.tsv"
    expression_metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{expression_metadata_file}"

    output_name = f"expression_data_{cancer}.pickle"
    params = {
        "expression_data_path": expression_data_path,
        "expression_metadata_location": expression_metadata_location,
        "output_name": output_name}
            
    expressionParser(params)

def performSNVanalysisInit():
    from tools import performSNVanalysis

    cancer = cfg.settings["TCGA-PROJECT"]
    username = cfg.settings["username"]
    metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
    metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"
    df_metadata = pd.read_csv(metadata_location, sep="\t")
    metadata_file_maf = f"TCGA-{cancer}-WXS-MAF-metadata.tsv"
    metadata_location_maf = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file_maf}"
    df_metadata_maf = pd.read_csv(metadata_location_maf, sep="\t")

    params = {
        "consolidated_results_path": cfg.settings["consolidated_results_path"],
        "df_metadata": df_metadata,
        "df_metadata_maf": df_metadata_maf,
        "MMBIR_THRESHOLD_LOW": cfg.settings["MMBIR_THRESHOLD_LOW"],
        "MMBIR_THRESHOLD_HIGH": cfg.settings["MMBIR_THRESHOLD_HIGH"],
        "min_concentration": cfg.settings["min_concentration"],
        "mmbir_fraction_high": cfg.settings["mmbir_fraction_high"],
        "mmbir_fraction_low": cfg.settings["mmbir_fraction_low"],
        "outputs_path": cfg.settings["outputs_path"]
        }

    performSNVanalysis(params)

def heatMapperInit():

    from tools import heatMapper

    file_path = sys.argv[2]

    df = pd.read_csv(file_path)

    #create a dictionary of positions for each interval from the dataframe.
    # the chr column is the key, and the values are all the positions in the sBirStart column for that chr
    positions = {}
    for chr in df["chr"].unique():
        positions[chr] = df[df["chr"] == chr]["iBirStart"].values

    #Rename chr24 to chrX, chr25 to chrY, and chr23 to chrM
    positions["chrM"] = positions.pop("chr23")
    positions["chrX"] = positions.pop("chr24")
    positions["chrY"] = positions.pop("chr25")

    #remove duplicate positions in the dictionary values (list of positions), say how many positions there are before and after removing duplicates
    for chr in positions:
        all_positions =  len(positions[chr])
        positions[chr] = list(dict.fromkeys(positions[chr]))
        deduplicated_positions = len(positions[chr])
        logging.info(f"Number of positions after removing duplicates for {chr}: {deduplicated_positions} (before: {all_positions}). Removed {all_positions - deduplicated_positions} duplicates.")

    heatMapper(positions, bandwidth=100000, tickspace=10000000, cmap="YlOrRd", save_path="outputs/heatmap.png")


if __name__ == "__main__":

    options_dict = {
        "help": help,
        "groupCases": groupCasesInit,
        "parseOutputs": parseOutputsInit,
        "createFullCancerTable": createFullCancerTableInit,
        "findCosmicGenes": findCosmicGenesInit,
        "performDiffExprAnalysis": performDiffExprAnalysisInit,
        "getMissingBams": getMissingBamsInit,
        "expressionParser": expressionParserInit,
        "performSNVanalysis": performSNVanalysisInit,
        "masked_snv_mv": masked_snv_mvInit,
        "heatMapper": heatMapperInit
        }

    if len(sys.argv) == 1:
        logging.error("No option selected. Use 'help' to see the list of available options.")
        sys.exit()

    option = sys.argv[1]

    if option in options_dict:
        options_dict[option]()

    else:
        logging.error(f"{option} is not a recognized option. Use 'help' to see the list of available options.")
