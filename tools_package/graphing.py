# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

from tools import *
import cancer_config as cfg
import os
import logging

logger = logging.getLogger(__name__)

# in LAML, must override default control_tumor=["Blood Derived Normal", "Primary Tumor"]
# to control_tumor = ["Solid Tissue Normal", "Primary Blood Derived Cancer - Peripheral Blood"]
# in the following functions:
#plot_blood_tumor_count_correlations_treshold_delta
#plot_blood_tumor_count_correlation
#plot_age_vs_count_correlation

def graphing(params, save=False, show=True):

    cancer = params["cancer"]
    metadata_location = params["metadata_location"]
    min_concentration = params["min_concentration"]
    staging = params["staging"]
    show = params["show_plots"]

    #load the consolidated results file, in future try migrating away from it to rely solely on mmbir_master_tables...
    df_consolidated=pd.read_csv(f"consolidated_results_{cancer}.tsv", sep="\t")
    df_metadata=pd.read_csv(metadata_location, sep="\t")

    df_consolidated = annotate_consolidated_results(df_consolidated, df_metadata)

    try:
        plot_differential_expression(cancer, save=save, show=show)
    except Exception as e:
        logging.error(f"Error in differential expression plot: {e}")

    try:
        plot_total_reads_vs_count(df_consolidated, count="Raw_Count", min_concentration=min_concentration, save=save, show=show)
    except Exception as e:
        logging.error(f"Error in total reads vs count plot: {e}")

    try:
        df_consolidated_trimmed = df_consolidated[df_consolidated["total_reads"] >= 100000000].copy()
        print("removing samples with fewer than 100000000 reads")
        plot_total_reads_vs_count(df_consolidated_trimmed, count="Raw_Count", min_concentration=min_concentration, save=save, show=show)
    except Exception as e:
        logging.error(f"Error in total reads vs count plot: {e}")

    try:
        df_wide = plot_blood_tumor_count_correlations_treshold_delta(df_consolidated, min_concentration=min_concentration, method="spearman", save=save, show=show) #method="spearman" or "pearson"
        plot_blood_tumor_count_correlation(df_wide, method="spearman", threshold=3000, save=save, show=show) #method="spearman" or "pearson"
    except Exception as e:
        logging.error(f"Error in blood tumor count correlation plot: {e}")

    try:
        plot_count_vs_concentration(df_consolidated, x_count="Raw_Count", save=save, show=show)
    except Exception as e:
        logging.error(f"Error in count vs concentration plot: {e}")

    filterset=["Blood Derived Normal","Primary Tumor"] #"Blood Derived Normal", "Primary Tumor"

    try:
        plot_concentration_raw_filtered(df_consolidated, filterset, hue="Concentration", save=save, show=show) #Concentration_bin
    except Exception as e:
        logging.error(f"Error in concentration raw filtered plot: {e}")

    try:
        plot_Sample_Type_counts(df_consolidated, filterset, min_concentration=min_concentration, cancer=cancer, save=save, show=show)
    except Exception as e:
        logging.error(f"Error in sample type counts plot: {e}")

    filterset=["Primary Tumor"] #"Blood Derived Normal", "Primary Tumor", control_tumor

    #staging='ajcc'/'figo'/'tumor_grade', adjust='early_late' or 'early_middle_late'
    try:
        plot_stage_vs_count(df_consolidated, filterset, staging=staging, x_count="Filtered_Count", min_concentration=min_concentration, adjust_staging='early_late', save=save, show=show)
    except Exception as e:
        logging.error(f"Error in stage vs count plot: {e}")
    try:
        plot_stage_vs_concentration(df_consolidated, filterset, staging=staging, x_count="Filtered_Count", min_concentration=min_concentration, adjust_staging='early_late', save=save, show=show)
    except Exception as e:
        logging.error(f"Error in stage vs concentration plot: {e}")

    if staging == 'ajcc':
        try:
            plot_ajcc_pathologic_n_vs_count(df_consolidated, filterset, x_count="Filtered_Count", min_concentration=min_concentration, save=save, show=show)
        except Exception as e:
            logging.error(f"Error in ajcc pathologic n vs count plot: {e}")
        try:
            plot_ajcc_pathologic_t_vs_count(df_consolidated, filterset, x_count="Filtered_Count", min_concentration=min_concentration, save=save, show=show)
        except Exception as e:
            logging.error(f"Error in ajcc pathologic t vs count plot: {e}")
    
    try:
        plot_age_vs_count_correlation(df_consolidated, count="Filtered_Count", min_concentration=min_concentration, method="spearman", save=save, show=show)
    except Exception as e:
        logging.error(f"Error in age vs count correlation plot: {e}")

    try:
        plot_age_vs_count_binned(df_consolidated, count="Filtered_Count", min_concentration=min_concentration, save=save, show=show)
    except Exception as e:
        logging.error(f"Error in age vs count binned plot: {e}")
    
    try:
        plot_age_vs_count_binned(df_consolidated, count="Raw_Count", min_concentration=min_concentration, save=save, show=show)
    except Exception as e:
        logging.error(f"Error in age vs count binned plot: {e}")

if __name__ == "__main__":

    cancer = cfg.settings["TCGA-PROJECT"]
    username = cfg.settings["username"]

    metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
    metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"

    params = {
        "cancer": cancer,
        "metadata_location": metadata_location,
        "min_concentration": 0,
        "staging": "ajcc",
        "save_plots": True,
        "show_plots": True,
    }

    #check if folder outputs/plots exists, if not, create it given that params["save"] is True
    if params["save_plots"]:
        if not os.path.exists("outputs/plots"):
            os.makedirs("outputs/plots")

    graphing(params, save=True)
