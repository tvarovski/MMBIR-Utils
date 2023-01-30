import pandas as pd
import os
import cancer_config as cfg

cancer = cfg.settings["TCGA-PROJECT"]
username = cfg.settings["username"]

metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"

df_metadata = pd.read_csv(metadata_location, sep="\t")
df_metadata = df_metadata.sort_values(by=["cases.0.case_id"])

#create a table for raw and filtered numbers of MMBIR events (headers and numbers)
output = os.system(f"echo 'sample	mmbir_events' > mmbir_counts_raw.tsv")
output = os.system(f"./count_events.sh outputs/raw >> mmbir_counts_raw.tsv")

output = os.system(f"echo 'sample	mmbir_events' > mmbir_counts_filtered.tsv")
output = os.system(f"./count_events.sh outputs/filtered >> mmbir_counts_filtered.tsv")

raw_df = pd.read_csv("mmbir_counts_raw.tsv", sep="\t")
raw_df["sample"]=raw_df["sample"].str[:-8]
filtered_df = pd.read_csv("mmbir_counts_filtered.tsv", sep="\t")
filtered_df["sample"]=filtered_df["sample"].str[:-29]

output = os.system(f"ls -d1 *-*-*-*-* > cases.txt")
print(f"cases exited with a code {output}")
cases_file = open("cases.txt", 'r')
cases = cases_file.readlines()
cases_file.close()


column_names = ["Case_ID", "Sample_Name", "Sample_Type", "Raw_Count", "Filtered_Count"]
output_df=pd.DataFrame(columns = column_names)

for caseID in cases:
  caseID=caseID.strip()
  case_files = df_metadata[df_metadata["cases.0.case_id"] == caseID]
  case_files = case_files[["file_name", "cases.0.samples.0.sample_type"]]
  output = os.system(f"ls --hide=*.* -1 {caseID} > {caseID}/ind_files.txt")

  samples_file = open(f"{caseID}/ind_files.txt", 'r')
  samples = samples_file.readlines()
  samples_file.close()
  
  parsing_list=[]
  for sample in samples:
    sample=f"{sample.strip()}.bam"
    df_sample = case_files[case_files["file_name"] == sample]
    sample_list = df_sample[["file_name", "cases.0.samples.0.sample_type"]].values.tolist()
    try:
      parsing_list.append(sample_list[0])
    except IndexError:
      print(f"Couldn't read data in sample list: {sample_list} of CaseID: {caseID}")
      exit()

  for i in parsing_list:

    raw_count_list = raw_df[raw_df["sample"]==i[0][:-4]]
    raw_count_list = raw_count_list["mmbir_events"].tolist()
    filtered_count_list = filtered_df[filtered_df["sample"]==i[0][:-4]]
    filtered_count_list = filtered_count_list["mmbir_events"].tolist()
    if len(filtered_count_list) == 0:
      filtered_count_list = 0
    append_me=pd.DataFrame({"Case_ID": caseID,
                            "Sample_Name": i[0],
                            "Sample_Type": i[1],
                            "Raw_Count": raw_count_list,
                            "Filtered_Count": filtered_count_list})

    output_df = pd.concat([append_me,output_df.loc[:]]).reset_index(drop=True)
  output = os.system(f"rm {caseID}/ind_files.txt")
output = os.system(f"rm cases.txt")

print(output_df)
output_df.to_csv("consolidated_results.tsv", sep="\t", index=False)
