import pandas as pd
import cancer_config as cfg
from tools import groupCases, parseOutputs

cancer = cfg.settings["TCGA-PROJECT"]
username = cfg.settings["username"]

metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"

df_metadata = pd.read_csv(metadata_location, sep="\t")

groupCases(df_metadata)

parseOutputs(df_metadata, consolidated_results_name="consolidated_results.tsv")
