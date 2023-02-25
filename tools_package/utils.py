import pandas as pd
import sys
import cancer_config as cfg
import logging


###DESCRIPTION###

'''

This is a script designed to be used from a command line which takes
a command line argument and runs a corresponding functionality from the
MMBIR-UTILS package

usage (from the command line):
    
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
        #!!!NOT WORKING!!!
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


'''

if __name__ == "__main__":

    #set logging level

    logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(name)s.%(funcName)s:%(message)s')
    logger = logging.getLogger(__name__)

    command = sys.argv[1]

    cancer = cfg.settings["TCGA-PROJECT"]
    username = cfg.settings["username"]

    metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
    metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"

    df_metadata = pd.read_csv(metadata_location, sep="\t")

    if command == "groupCases":

        from tools import groupCases
        groupCases(df_metadata)

    elif command == "parseOutputs":

        from tools import parseOutputs
        parseOutputs(df_metadata, consolidated_results_name="consolidated_results.tsv")

    elif command == "createFullCancerTable":

        from tools import createFullCancerTable
        manifest_file = f"TCGA-{cancer}-WXS-BAM-manifest.tsv"

        params={
            "df_metadata": df_metadata,
            "manifest_location": f"/Users/{username}/MMBIR_Databases/TCGA/{manifest_file}",
            "missing_manifest_output_name": f"missing_bams_{cancer}_manifest.tsv",
            "log": False,
            "output_raw_name": "raw_mmbir_results_master_df.csv",
            "output_filtered_name": "filtered_mmbir_results_master_df.csv"
            }
            
        createFullCancerTable(params)

    elif command == "findCosmicGenes":

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

    elif command == "performDiffExprAnalysis":

        from tools import performDiffExprAnalysis
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
        
    elif command == "getMissingBams":

        from tools import getMissingBams
        manifest_file = f"TCGA-{cancer}-WXS-BAM-manifest.tsv"
        manifest_location = f"/Users/{username}/MMBIR_Databases/TCGA/{manifest_file}"
        missing_manifest_output_name = f"missing_bams_{cancer}_manifest.tsv"

        getMissingBams(df_metadata, manifest_location, missing_manifest_output_name)

    elif command == "expressionParser":
        
        from tools import expressionParser

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

    elif command == "masked_snv_mv":

        from tools import masked_snv_mv

        #need to add a sys.argv[2] for the snvs_loc
        snvs_loc = sys.argv[2]

        metadata_file_maf = f"TCGA-{cancer}-WXS-MAF-metadata.tsv"
        metadata_location_maf = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file_maf}"
        df_metadata = pd.read_csv(metadata_location_maf, sep="\t")

        masked_snv_mv(df_metadata, snvs_loc=snvs_loc)

    elif command == "performSNVanalysis":

        from tools import performSNVanalysis

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

    else:
        logging.error(f"{command} is not a recognized option")

