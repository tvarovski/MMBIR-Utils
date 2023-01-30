import pandas as pd
import os
import cancer_config as cfg

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

def get_case_bam_files(case_id):
    #this function returns a list of bam files for a given case_id
    df_case = df_metadata[df_metadata["cases.0.case_id"]==case_id]
    bam_files = df_case["file_name"].tolist()
    return bam_files

def get_bam_file_metadata(bam_file):
    #this function returns metadata associated with a given bam file
    df_bam = df_metadata[df_metadata["file_name"]==bam_file]
    #return a dictionary of metadata
    return df_bam

def check_for_missing_bams():
    #go through all the cases and create a list of all the files
    cases_list = []
    for index, row in df_metadata.iterrows():
        case_id = row['cases.0.case_id']
        cases_list.append(case_id)

    #remove duplicates
    cases_list = list(dict.fromkeys(cases_list))

    #check if the file exists as a directory
    found_cases=0
    found_bams=0
    missing_bams = []
    missing_cases = []
    for case_id in cases_list:
        #get current directory
        current_dir = os.getcwd()
        case_dir = f"{current_dir}/{case_id}"
        if not os.path.exists(case_dir):
            print(f"WARNING: The directory: {case_dir} does not exist")
            missing_cases.append(case_dir)
        else:
            found_cases+=1
            #enter the directory
            os.chdir(case_dir)
            #check what directories are in the directory
            dir_list = os.listdir()
            
            #get the bam files for the case
            case_bam_files = get_case_bam_files(case_id)
            #check if the bam files are in the directory
            for bam_file in case_bam_files:
                if bam_file.strip(".bam") not in dir_list:
                    print(f"WARNING: The bam file: {bam_file} does not exist in: {case_dir}")
                    missing_bams.append(bam_file)
                else:
                    found_bams+=1
            #go back to the current directory
            os.chdir(current_dir)
    print(f"Found {found_cases} cases out of {len(cases_list)}")
    print(f"Found {found_bams} bam files out of {len(df_metadata)}")
    print("\nMissing bams:\n")
    print(missing_bams)
    print("Missing cases:\n")
    print(missing_cases)
    return([missing_bams, missing_cases])

def create_missing_bams_manifest(missing_files, manifest_location, missing_manifest_output_name):

    if len(missing_files[0]) == 0 & len(missing_files[0]) == 0:
        print("All good! All files present!")
        return
    else:
        # create an empty, placeholder dataframe
        df_out = pd.DataFrame()
        #open the manifest file
        df = pd.read_csv(manifest_location, sep="\t")
        if len(missing_files[0]) > 0:
            df_missing_bams = df[df["filename"].isin(missing_files[0])]
            df_out = pd.concat([df_out, df_missing_bams], axis=0)
            
        if len(missing_files[1]) > 0:
            #open the metadata file to retreive filenames associated with missing caseIDs
            df_metadata = pd.read_csv(metadata_location, sep="\t")
            df_metadata_missing_cases = df_metadata[df_metadata["cases.0.case_id"].isin(missing_files[1])]

            #get the manifest rows containing only files from missing cases
            df_missing_cases = df[df["filename"].isin(df_metadata_missing_cases["file_name"])]
            df_out = pd.concat([df_out, df_missing_cases], axis=0)
        
        #save the df_out to the new output_file
        df_out.to_csv(missing_manifest_output_name, sep="\t", index=False)
        return(df_out)

def create_mmbir_results_master_df(filtered=False, log=False):
    # create a master dataframe for raw mmbir results

    results_master_df = pd.DataFrame()
    if filtered:
        result_dir="filtered"
    else:
        result_dir="raw"
    #get current directory
    current_dir = os.getcwd()
    os.chdir(f"{current_dir}/outputs/{result_dir}")
    #list the files in the directory
    dir_list = os.listdir()
    #go through the files in the directory
    for file_name in dir_list:

        file_name = file_name.split("/")[-1]
        bam_name = file_name.replace("_raw.txt", "")
        bam_name = bam_name.replace("_allchrmmbir_filtered_all.txt", "")
        bam_name += ".bam"

        bam_metadata = get_bam_file_metadata(bam_name)
        #print(bam_metadata)

        if len(bam_metadata) == 0:
            print(f"WARNING: No metadata found for: {bam_name}, skipping")
            continue
        elif len(bam_metadata) > 1:
            print(f"WARNING: Multiple metadata found for: {bam_name}, skipping")
            continue

        #get column names of the bam_metadata dataframe
        column_names = bam_metadata.columns.tolist()

        #read in the filtered mmbir results
        filtered_mmbir_results_df = pd.read_csv(file_name, sep="\t")

        #add the bam_metadata to the filtered_mmbir_results_df, rename the column names for brevity
        for column_name in column_names:
            new_column_name = column_name.replace("cases.0.", "")
            new_column_name = new_column_name.replace("demographic.", "")
            new_column_name = new_column_name.replace("diagnoses.0.", "")
            new_column_name = new_column_name.replace("project.", "")
            new_column_name = new_column_name.replace("samples.0.", "")
            new_column_name = new_column_name.replace("portions.0.", "")
            new_column_name = new_column_name.replace("analytes.0.", "")

            #add the bam metadata to the dataframe
            filtered_mmbir_results_df[new_column_name] = bam_metadata[column_name].tolist()[0]
        
        #concat the filtered mmbir results to the master dataframe
        results_master_df = pd.concat([results_master_df, filtered_mmbir_results_df], axis=0)
        if log:
            print(f"Finished processing: {bam_name}")
    
    #go back to the current directory
    os.chdir(current_dir)
    return results_master_df

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

