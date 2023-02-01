import pandas as pd
import cancer_config as cfg
from tools import createFullCancerTable


if __name__ == "__main__":

    cancer = cfg.settings["TCGA-PROJECT"]
    username = cfg.settings["username"]

    metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
    manifest_file = f"TCGA-{cancer}-WXS-BAM-manifest.tsv"
    metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"
    df_metadata = pd.read_csv(metadata_location, sep="\t")

    params={
        "df_metadata": df_metadata,
        "manifest_location": f"/Users/{username}/MMBIR_Databases/TCGA/{manifest_file}",
        "missing_manifest_output_name": f"missing_bams_{cancer}_manifest.tsv",
        "log": False,
        "output_raw_name": "raw_mmbir_results_master_df.csv",
        "output_filtered_name": "filtered_mmbir_results_master_df.csv"
        }

    createFullCancerTable(params)
