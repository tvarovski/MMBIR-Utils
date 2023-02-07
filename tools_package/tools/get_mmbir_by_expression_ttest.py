import pandas as pd
import cancer_config as cfg
from tools import performDiffExprAnalysis

if __name__ == "__main__":
    
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
        "min_concentration": cfg.settings["min_concentration"]
    }


    performDiffExprAnalysis(params)

