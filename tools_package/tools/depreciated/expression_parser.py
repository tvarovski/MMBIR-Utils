import cancer_config as cfg
from tools import expressionParser


if __name__ == "__main__":

    # set environment variables
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
        "output_name": output_name
    }

    expressionParser(params)

