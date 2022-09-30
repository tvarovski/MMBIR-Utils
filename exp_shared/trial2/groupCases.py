import shutil 
import os
import pandas as pd


metadata_file_path="/Users/twarowski/MMBIR_Databases/TCGA/TCGA-OV-WXS-BAM-metadata.tsv"

df_metadata = pd.read_csv(metadata_file_path, sep="\t")


mypath=os.path.abspath(os.getcwd())

onlydirs = [d for d in os.listdir(mypath) if os.path.isdir(os.path.join(mypath, d))]

print(f"Found {len(onlydirs)} directories: {onlydirs}")



for current_name in onlydirs:

  current_file_entry = df_metadata[df_metadata["file_name"] == f"{current_name}.bam"]
  current_case_id = current_file_entry["cases.0.case_id"].tolist()
  if len(current_case_id) == 0:
    continue
  if len(current_case_id) != 1:
    print("Retreived unexpected number of case_id for sample metadata. Exitting...")
    print(f"Case IDs: {current_case_id}")
    exit()
  current_case_id_samples = df_metadata[df_metadata["cases.0.case_id"] == current_case_id[0]]

  current_case_id_samples_list = current_case_id_samples["file_name"].tolist()

  for sample in current_case_id_samples_list:
    current_case_id_path = os.path.join(mypath, current_case_id[0])
    if not os.path.exists(current_case_id_path):
      os.makedirs(current_case_id_path)

    if os.path.exists(os.path.join(mypath, sample[:-4])):
      shutil.move(sample[:-4], current_case_id_path)
      print(f"Moved {sample[:-4]} to {current_case_id_path}")
