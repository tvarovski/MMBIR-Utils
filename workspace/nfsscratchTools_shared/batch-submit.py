import os

manifest_path="/Users/twarowski/MMBIR_Databases/TCGA/TCGA-BRCA-WXS-BAM-manifest.tsv"
slice_path="/nfsscratch/twarowski/TCGA-BRCA/slice-0/"


df_manifest = pd.read_csv(manifest_path, sep="\t")
job_params_list = df_manifest[["id", "filename"]].tolist()

for id, filename in job_params_list:
  filename = filename.strip(".bam")
  sample_path=slice_path+id
  os.system(f"qsub smk_pipeline.job {sample_path} {filename}")
