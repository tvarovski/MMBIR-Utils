# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

import os
import pandas as pd
import subprocess

SLICE_NUM=0
USER="twarowski"

manifest_path=os.path.join("/Users", USER, "MMBIR_Databases", "TCGA", "manifest_slices", f"TCGA-COAD-WXS-BAM-manifest-slice-{SLICE_NUM}.txt")
slice_path=os.path.join("/nfsscratch", USER, "TCGA-COAD", f"slice-{SLICE_NUM}")

df_manifest = pd.read_csv(manifest_path, sep="\t")
job_params_list = df_manifest[["id", "filename"]].values.tolist()

for id, filename in job_params_list:
    filename = filename[:-4]
    sample_path=os.path.join(slice_path, id)
    command = f"qsub smk_pipeline.job {sample_path} {filename}"
    subprocess.run(command, shell=True)



'''
import os
import pandas as pd


SLICE_NUM=0
USER="twarowski"

manifest_path=f"/Users/twarowski/MMBIR_Databases/TCGA/manifest_slices/TCGA-COAD-WXS-BAM-manifest-slice-{SLICE_NUM}.txt"
slice_path=f"/nfsscratch/{USER}/TCGA-COAD/slice-{SLICE_NUM}/"


df_manifest = pd.read_csv(manifest_path, sep="\t")
job_params_list = df_manifest[["id", "filename"]].values.tolist()

for id, filename in job_params_list:
    filename = filename[:-4]
    sample_path=slice_path+id
    os.system(f"qsub smk_pipeline.job {sample_path} {filename}")
'''