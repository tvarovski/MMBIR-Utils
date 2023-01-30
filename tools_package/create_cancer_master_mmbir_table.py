import pandas as pd
import os
import cancer_config as cfg
from tools import get_case_bam_files, get_bam_file_metadata, check_for_missing_bams, create_missing_bams_manifest, create_mmbir_results_master_df

cancer = cfg.settings["TCGA-PROJECT"]
username = cfg.settings["username"]

metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"

manifest_file = f"TCGA-{cancer}-WXS-BAM-manifest.tsv"
manifest_location = f"/Users/{username}/MMBIR_Databases/TCGA/{manifest_file}"

output_raw_name = "raw_mmbir_results_master_df.csv"
output_filtered_name = "filtered_mmbir_results_master_df.csv"

missing_manifest_output_name = f"missing_bams_{cancer}_manifest.tsv"

df_metadata = pd.read_csv(metadata_location, sep="\t")


# create a main function
def main(log=False):
    my_dir=os.getcwd()
    #os.chdir(f"{my_dir}/finished")
    #check for missing bams
    missing_files = check_for_missing_bams()
    #pass the list, if missing create manifest for redownload
    create_missing_bams_manifest(missing_files, manifest_location, missing_manifest_output_name)
    #create the master dataframe for the raw mmbir results
    raw_mmbir_results_master_df = create_mmbir_results_master_df(filtered=False)
    print(f"Finished creating the master dataframe for the raw mmbir results")
    #create the master dataframe for the filtered mmbir results
    filtered_mmbir_results_master_df = create_mmbir_results_master_df(filtered=True)
    print(f"Finished creating the master dataframe for the filtered mmbir results")

    os.chdir(my_dir)
    #save the master dataframes to csv
    raw_mmbir_results_master_df.to_csv(output_raw_name, index=False)
    filtered_mmbir_results_master_df.to_csv(output_filtered_name, index=False)
    print(f"Finished saving the master dataframes to csv")

main()

