from tools import *
import cancer_config as cfg
import os


@fancy_status
def plot_differential_expression(cancer, save=False, show=True):
    import numpy as np
    path = f"outputs/ttest_results_{cancer}_minconc0_bh_corrected.tsv"

    df = pd.read_csv(path, sep="\t")

    #transform the p-value to -log10
    df["-log10(p-value)"] = -np.log10(df["p-value"])
    
    #transform the fold change to log2
    df["log2(fold change)"] = np.log2(df["fold-change"])

    #if log2(fold change) is greater than 10, set it to 10
    #if flog2(fold change) is less than -10, set it to -10
    df["log2(fold change)"] = df["log2(fold change)"].apply(lambda x: 6 if x > 6 else x)
    df["log2(fold change)"] = df["log2(fold change)"].apply(lambda x: -6 if x < -6 else x)

    #plot the -log10(p-value) vs the log2(fold change), color the points blue if the p-value is less than 0.05
    #and red if the p-value is greater than 0.05
    #also, add a line at -log10(p-value) = 1.3, which is the threshold for significance
    df["significant"] = df["p-value"] < 0.05

    sns.scatterplot(x="log2(fold change)", y="-log10(p-value)", data=df, hue="significant")
    plt.axhline(y=1.3, color="red", linestyle="--")
    plt.axvline(x=0, color="black", linestyle="--")

    #rename the x-axis and y-axis
    plt.xlabel("log2(fold change)", fontsize=20, fontweight="bold")
    plt.ylabel("-log10(p-value)", fontsize=20, fontweight="bold")

    #add the title
    plt.title(f"Differential expression of {cancer}", fontsize=20, fontweight="bold")
    plt.tight_layout()
    
    #add the a vertical line at log2(fold change) = 0, which is the threshold for differential expression

    if save:
        plt.savefig(f"outputs/plots/{cancer}_differential_expression.png", dpi=600)

    if show:
        plt.show()
    else:
        plt.close()

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

    plot_differential_expression(cancer, save=save, show=show)

    plot_total_reads_vs_count(df_consolidated, count="Raw_Count", min_concentration=min_concentration, save=save, show=show)

    df_consolidated = df_consolidated[df_consolidated["total_reads"] >= 100000000]
    print("removing samples with fewer than 100000000 reads")

    plot_total_reads_vs_count(df_consolidated, count="Raw_Count", min_concentration=min_concentration, save=save, show=show)

    df_wide = plot_blood_tumor_count_correlations_treshold_delta(df_consolidated, min_concentration=min_concentration, method="spearman", save=save, show=show) #method="spearman" or "pearson"
    plot_blood_tumor_count_correlation(df_wide, method="spearman", threshold=3000, save=save, show=show) #method="spearman" or "pearson"

    plot_count_vs_concentration(df_consolidated, x_count="Raw_Count", save=save, show=show)

    filterset=["Blood Derived Normal","Primary Tumor"] #"Blood Derived Normal", "Primary Tumor"

    plot_concentration_raw_filtered(df_consolidated, filterset, hue="Concentration", save=save, show=show) #Concentration_bin
    plot_Sample_Type_counts(df_consolidated, filterset, min_concentration=min_concentration, cancer=cancer, save=save, show=show)

    filterset=["Primary Tumor"] #"Blood Derived Normal", "Primary Tumor"

    #staging='ajcc'/'figo'/'tumor_grade', adjust='early_late' or 'early_middle_late'
    plot_stage_vs_count(df_consolidated, filterset, staging=staging, x_count="Filtered_Count", min_concentration=min_concentration, adjust_staging='early_late', save=save, show=show)
    plot_stage_vs_concentration(df_consolidated, filterset, staging=staging, x_count="Filtered_Count", min_concentration=min_concentration, adjust_staging='early_late', save=save, show=show)

    if staging == 'ajcc':
        plot_ajcc_pathologic_n_vs_count(df_consolidated, filterset, x_count="Filtered_Count", min_concentration=min_concentration, save=save, show=show)
        plot_ajcc_pathologic_t_vs_count(df_consolidated, filterset, x_count="Filtered_Count", min_concentration=min_concentration, save=save, show=show)

    plot_age_vs_count_correlation(df_consolidated, count="Filtered_Count", min_concentration=min_concentration, method="spearman", save=save, show=show)
    
    plot_age_vs_count_binned(df_consolidated, count="Filtered_Count", min_concentration=min_concentration, save=save, show=show)
    plot_age_vs_count_binned(df_consolidated, count="Raw_Count", min_concentration=min_concentration, save=save, show=show)

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
