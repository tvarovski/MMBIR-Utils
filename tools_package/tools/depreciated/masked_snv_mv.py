import pandas as pd
import cancer_config as cfg
from tools import masked_snv_mv

cancer = cfg.settings["TCGA-PROJECT"]
username = cfg.settings["username"]

snvs_metadata=f"/Users/{username}/MMBIR_Databases/TCGA/TCGA-{cancer}-WXS-MAF-metadata.tsv"
snvs_loc=f"/nfsscratch/{username}/TCGA-{cancer}/Masked_MAF"

df_metadata = pd.read_csv(snvs_metadata, sep="\t")

masked_snv_mv(df_metadata)