import pandas as pd
import sys
import cancer_config as cfg
from tools import groupCases, parseOutputs, masked_snv_mv

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

else:
    print(f"{command} is not a recognized option")