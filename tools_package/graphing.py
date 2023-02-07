from tools import *
import cancer_config as cfg
from tools import graphing


def graphing(params):

    cancer = params["cancer"]
    metadata_location = params["metadata_location"]
    min_concentration = params["min_concentration"]

    #load the consolidated results file, in future try migrating away from it to rely solely on mmbir_master_tables...
    df_consolidated=pd.read_csv(f"consolidated_results_{cancer}.tsv", sep="\t")
    df_metadata=pd.read_csv(metadata_location, sep="\t")

    df_consolidated = annotate_consolidated_results(df_consolidated, df_metadata)

    df_wide = plot_blood_tumor_count_correlations_treshold_delta(df_consolidated, min_concentration=min_concentration, method="spearman") #method="spearman" or "pearson"
    plot_blood_tumor_count_correlation(df_wide, method="spearman", threshold=3000)

    plot_count_vs_concentration(df_consolidated, x_count="Raw_Count")

    filterset=["Blood Derived Normal","Primary Tumor"] #"Blood Derived Normal", "Primary Tumor"

    plot_concentration_raw_filtered(df_consolidated, filterset, hue="Concentration") #Concentration_bin
    plot_Sample_Type_counts(df_consolidated, filterset, min_concentration=min_concentration, cancer=cancer)

    filterset=["Primary Tumor"] #"Blood Derived Normal", "Primary Tumor"

    #staging='ajcc' or 'figo', adjust='early_late' or 'early_middle_late'
    plot_stage_vs_count(df_consolidated, filterset, staging='ajcc', x_count="Filtered_Count", min_concentration=min_concentration, adjust_staging='early_late')
    plot_stage_vs_concentration(df_consolidated, filterset, staging='ajcc', x_count="Filtered_Count", min_concentration=min_concentration, adjust_staging='early_late')

    plot_ajcc_pathologic_n_vs_count(df_consolidated, filterset, x_count="Filtered_Count", min_concentration=min_concentration)
    plot_ajcc_pathologic_t_vs_count(df_consolidated, filterset, x_count="Filtered_Count", min_concentration=min_concentration)

    plot_age_vs_count_correlation(df_consolidated, count="Filtered_Count", min_concentration=min_concentration, method="spearman")
    
    plot_age_vs_count_binned(df_consolidated, count="Filtered_Count", min_concentration=min_concentration)
    plot_age_vs_count_binned(df_consolidated, count="Raw_Count", min_concentration=min_concentration)


if __name__ == "__main__":

    cancer = cfg.settings["TCGA-PROJECT"]
    username = cfg.settings["username"]

    metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
    metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"

    params = {
        "cancer": cancer,
        "metadata_location": metadata_location,
        "min_concentration": 0
    }

    graphing(params)
