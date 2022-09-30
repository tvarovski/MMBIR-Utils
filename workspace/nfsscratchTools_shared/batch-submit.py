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
