import pandas as pd
import os
import cancer_config as cfg

cancer = cfg.settings["TCGA-PROJECT"]
username = cfg.settings["username"]

snvs_metadata=f"/Users/{username}/MMBIR_Databases/TCGA/TCGA-{cancer}-WXS-MAF-metadata.tsv"
snvs_loc=f"/nfsscratch/{username}/TCGA-{cancer}/Masked_MAF"

df_metadata = pd.read_csv(snvs_metadata, sep="\t")

output = os.system(f"ls -d1 *-*-*-*-* > cases.txt")
print(f"cases exited with a code {output}")
cases_file = open("cases.txt", 'r')
cases = cases_file.readlines()
cases_file.close()


for caseID in cases:
    
  caseID=caseID.strip()
  case_snvs = df_metadata[df_metadata["cases.0.case_id"] == caseID]
  case_files = case_snvs[["file_name", "id"]].values.tolist()
  for file_name, id in case_files:
    output = os.system(f"cp {snvs_loc}/{id}/{file_name} {caseID}")
    if output == 1:
      print(f"there is a problem with {snvs_loc}/{id}/{file_name} {caseID}")
    output = os.system(f"yes n | gunzip {caseID}/{file_name}")
