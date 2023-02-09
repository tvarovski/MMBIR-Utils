import pandas as pd
import cancer_config as cfg
from tools import performSNVanalysis


if __name__ == "__main__":
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
        "min_concentration": 0
    }


    performSNVanalysis(params)