import pandas as pd
import sys
import cancer_config as cfg
from tools import groupCases, parseOutputs, masked_snv_mv, createFullCancerTable, findCosmicGenes

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
        #filteres out and prints out genes in the Cosmic_DB that
        #are present in the specified DF path (sys.argv[2])

'''


command = sys.argv[1]

cancer = cfg.settings["TCGA-PROJECT"]
username = cfg.settings["username"]

metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"

df_metadata = pd.read_csv(metadata_location, sep="\t")

if command == "groupCases":
    groupCases(df_metadata)

elif command == "parseOutputs":
    parseOutputs(df_metadata, consolidated_results_name="consolidated_results.tsv")

elif command == "masked_snv_mv":
    masked_snv_mv(df_metadata)

elif command == "createFullCancerTable":

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
    

else:
    print(f"{command} is not a recognized option")