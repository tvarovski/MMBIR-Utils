from tools import *
import cancer_config as cfg


if __name__ == "__main__":

    cancer = cfg.settings["TCGA-PROJECT"]
    username = cfg.settings["username"]

    metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
    metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"

    #load the consolidated results file
    df_consolidated=pd.read_csv(f"consolidated_results_{cancer}.tsv", sep="\t")
    df_metadata=pd.read_csv(metadata_location, sep="\t")

    df_consolidated = annotate_consolidated_results(df_consolidated, df_metadata)

    df_wide = plot_blood_tumor_count_correlations_treshold_delta(df_consolidated, min_concentration=0.5, method="spearman") #method="spearman" or "pearson"
    plot_blood_tumor_count_correlation(df_wide, method="spearman", threshold=3000)

    plot_count_vs_concentration(df_consolidated, x_count="Raw_Count")

    filterset=["Blood Derived Normal","Primary Tumor"] #"Blood Derived Normal", "Primary Tumor"

    plot_concentration_raw_filtered(df_consolidated, filterset, hue="Concentration") #Concentration_bin
    plot_Sample_Type_counts(df_consolidated, filterset, min_concentration=0.5, cancer=cancer)

    filterset=["Primary Tumor"] #"Blood Derived Normal", "Primary Tumor"

    #staging='ajcc' or 'figo', adjust='early_late' or 'early_middle_late'
    plot_stage_vs_count(df_consolidated, filterset, staging='figo', x_count="Filtered_Count", min_concentration=0.5, adjust_staging='early_late')
    plot_stage_vs_concentration(df_consolidated, filterset, staging='figo', x_count="Filtered_Count", min_concentration=0.5, adjust_staging='early_late')

    plot_ajcc_pathologic_n_vs_count(df_consolidated, filterset, x_count="Filtered_Count", min_concentration=0.5)
    plot_ajcc_pathologic_t_vs_count(df_consolidated, filterset, x_count="Filtered_Count", min_concentration=0.5)

    plot_count_vs_age_correlation(df_consolidated, count="Filtered_Count", min_concentration=0.5, method="spearman")