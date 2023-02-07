# Plot number of MMBIR events in Tumor vs number of MMBIR events in Control

from importlib.metadata import distribution
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from tools import fancy_status

@fancy_status
def annotate_consolidated_results(df_consolidated, df_metadata):
    '''function annotate_consolidated_results() takes two dataframes (df_consolidated and df_metadata)
    as input, performs a join on specific columns, performs some data processing, and returns the updated dataframe.'''

    #join df_metadata and df_consolidated on the file_name column for df_metadata and for Sample_Name for df_consolidated
    df_consolidated=pd.merge(df_consolidated, df_metadata, left_on="Sample_Name", right_on="file_name")

    df_consolidated["age_at_collection"] = df_consolidated["cases.0.diagnoses.0.age_at_diagnosis"] + df_consolidated["cases.0.samples.0.days_to_collection"]
    #divide age_at_collection column by 365 to get the age in years
    df_consolidated["age_at_collection"] = df_consolidated["age_at_collection"]/365

    df_consolidated["age_at_collection"] = df_consolidated["cases.0.demographic.days_to_birth"]
    #divide age_at_collection column by 365 to get the age in years
    df_consolidated["age_at_collection"] = df_consolidated["age_at_collection"]/365 * -1

    #rename the cases.0.samples.0.portions.0.analytes.0.aliquots.0.concentration to concentration
    df_consolidated.rename(columns={"cases.0.samples.0.portions.0.analytes.0.aliquots.0.concentration": "Concentration"}, inplace=True)

    return df_consolidated

@fancy_status
def plot_blood_tumor_count_correlations_treshold_delta(df_consolidated, min_concentration=0.5, method="spearman"):
    '''function plot_blood_tumor_count_correlations_treshold_delta() takes the output from
    annotate_consolidated_results() function as input and plots the correlations between 
    blood-derived normal and primary tumor sample counts using different thresholds. 
    The correlations are calculated using either the spearman or pearson method. 
    The plot shows the R-squared values for different thresholds and draws a horizontal line at 0.7.'''

    df_consolidated = df_consolidated[df_consolidated["Concentration"]>=min_concentration]

    #remove the Sample_Name column
    df_consolidated.drop(columns=["Sample_Name"], inplace=True)

    #groupby the Case_ID and Sample_Type columns and get the max of the Raw_Count and Filtered_Count columns
    agg_dict={"Raw_Count": ['max'],
                "Filtered_Count": ['max']}
    df_agg = df_consolidated.groupby(["Case_ID", "Sample_Type"]).agg(agg_dict).reset_index()

    # flatten the columns
    df_agg.columns = ['_'.join(col).strip() for col in df_agg.columns.values]

    # get Sample_Type that is only "Blood Derived Normal" or "Primary Tumor"
    df_agg = df_agg[df_agg["Sample_Type_"].isin(["Blood Derived Normal", "Primary Tumor"])]

    # transform the df to the wide format
    df_wide = df_agg.pivot(index="Case_ID_", columns="Sample_Type_", values="Filtered_Count_max").reset_index() #Filtered_Count_max

    # rename the columns
    df_wide.columns = ["Case_ID", "Blood_Derived_Normal", "Primary_Tumor"]

    #remove NAN rows
    df_wide = df_wide.dropna()
    print(df_wide)
    
    # find correlation for range of thresholds from 0 to 3000 with the step of 50
    thresholds = range(0,3000, 50)

    correlations = []

    for threshold in thresholds:
        df_wide_thr = df_wide[(df_wide["Blood_Derived_Normal"] <= threshold) & (df_wide["Primary_Tumor"] <= threshold)]
        r_squared = df_wide_thr.corr(method=method)["Blood_Derived_Normal"]["Primary_Tumor"]
        correlations.append([threshold, r_squared])

    for correlation in correlations:
        print(f"Threshold: {correlation[0]}, R-squared: {correlation[1]}")

    # convert to pandas dataframe
    df_correlations = pd.DataFrame(correlations, columns=["Threshold", "R-squared"])

    # set the plot context
    sns.set_context("talk")

    # plot the scatter plot of the Threshold and R-squared
    sns.lineplot(x="Threshold", y="R-squared", data=df_correlations)

    #draw a horizontal line at 0.7
    plt.axhline(y=0.7, color='r', linestyle=':')

    # rename the x-axis and y-axis
    plt.xlabel("Max Patient MMBIR cut-off")
    plt.ylabel("Correlation (R-squared)")

    plt.show()

    return df_wide

@fancy_status
def plot_blood_tumor_count_correlation(df_wide, method="spearman", threshold=3000):
    '''function plot_blood_tumor_count_correlation() plots the correlations between
    blood-derived normal and primary tumor sample counts using a threshold.
    The correlations are calculated using either the spearman or pearson method.
    The plot shows the R-squared values for different thresholds and draws a horizontal line at 0.7.'''

    df_wide_size_initial=len(df_wide)
    df_wide = df_wide[(df_wide["Blood_Derived_Normal"] <= threshold) & (df_wide["Primary_Tumor"] <= threshold)]
    df_wide_size_threshold=len(df_wide)

    print(f"Initial size: {df_wide_size_initial}, Threshold size: {df_wide_size_threshold}, percentage: {df_wide_size_threshold/df_wide_size_initial*100}")

    # plot the Blood_Derived_Normal vs Primary_Tumor on a scatter plot and find R-squared
    sns.regplot(x="Primary_Tumor",y="Blood_Derived_Normal", data=df_wide)
    r_squared = df_wide.corr(method=method)["Blood_Derived_Normal"]["Primary_Tumor"]
    # add the R-squared to the plot title, round the R-squared to 2 decimal places
    plt.title(f"R-squared: {r_squared:.2f}, Max Patient MMBIR cut-off: {threshold}")

    #rename the x-axis
    plt.xlabel("Primary Tumor MMBIR count")
    #rename the y-axis
    plt.ylabel("Blood Derived Normal MMBIR count")
    # calculate the R-squared
    print(f"R-squared is: {r_squared}")

    plt.show()

@fancy_status
def plot_count_vs_concentration(df_consolidated, x_count="Raw_Count"):
    '''function plot_count_vs_concentration() plots the concentration vs the MMBIR count.
    The plot shows the MMBIR count on the x-axis and the concentration on the y-axis.'''

    # plot the number of MMBIR events in of "Raw_Count" and "cases.0.samples.0.portions.0.analytes.0.aliquots.0.concentration"
    sns.scatterplot(x=x_count,y="Concentration", data=df_consolidated, alpha=0.5, hue="Sample_Type")
    # set log scale for x-axis
    plt.xscale("log")
    #rename the x-axis
    plt.xlabel("Raw MMBIR Signature Count")
    #rename the y-axis
    plt.ylabel("Aliquot Concentration")
    plt.show()

@fancy_status
def plot_concentration_raw_filtered(df_consolidated, filterset, hue="Concentration"):
    '''function plot_concentration_raw_filtered() plots the concentration vs the MMBIR count.
    The plot shows the Filtered MMBIR count on the x-axis and the Raw MMBIR count on y-axis.'''


    df_consolidated = df_consolidated[df_consolidated["Sample_Type"].isin(filterset)]

    print("Mean and standard deviation of the Filtered_Count:")
    print(df_consolidated.groupby("Sample_Type").mean()["Filtered_Count"])
    print(df_consolidated.groupby("Sample_Type").std()["Filtered_Count"])

    print("Cases with concentration of 0.5:")
    print(df_consolidated[df_consolidated["Concentration"]==.5])

    if hue == "Concentration_bin":
        #bin the concentration into 2 groups for plotting 2 discrete colors based on the concentration
        df_consolidated["Concentration_bin"] = "none"
        df_consolidated["Concentration_bin"] = df_consolidated["Concentration"].apply(lambda x: "low" if x < .5 else "high")

    # set the plot context
    sns.set_context("poster")
    sns.scatterplot(x="Filtered_Count", y="Raw_Count", data=df_consolidated, alpha=0.5, hue=hue, palette="flare")

    # set both x- and y-axis to log scale
    plt.xscale("log")
    plt.yscale("log")

    #rename the x-axis and y-axis
    plt.xlabel("Filtered MMBIR Signature Count")
    plt.ylabel("Raw MMBIR Signature Count")

    plt.show()

@fancy_status
def plot_Sample_Type_counts(df_consolidated, filterset, min_concentration=0.5, cancer="cancer"):
    '''function plot_Sample_Type_counts() plots the concentration vs the MMBIR count.
    The plot shows the MMBIR count on the x-axis and the concentration on the y-axis.'''

    df_figure = df_consolidated[df_consolidated["Sample_Type"].isin(filterset)]
    df_figure = df_figure[df_figure["Concentration"] >= min_concentration]
    
    sns.set_context("poster")
    #set color for each sample type
    color_dict = {"Primary Tumor": "tab:blue", "Blood Derived Normal": "tab:orange"}
    sns.histplot(data=df_figure, x="Raw_Count", bins=400, multiple="stack", hue="Sample_Type", palette=color_dict) #multiple="stack", hue="Sample_Type",
    
    #make labels bold
    plt.xlabel("MMBIR Signature Count", fontweight="bold")
    plt.ylabel("Number of tumor samples", fontweight="bold")
    
    #add title
    #show the length of the df_figure on the title
    length=len(df_figure)
    plt.title(f"{cancer} Cancer (N={length})", fontweight="bold")

    #don't show the legend
    #plt.legend().remove()

    #set x-axis to go to 1080
    #plt.xlim(0,1050)
    
    #make tick labels bold
    plt.xticks(fontweight="bold", fontsize=20)
    plt.yticks(fontweight="bold", fontsize=20)
    plt.show()

@fancy_status
def plot_stage_vs_count(df_consolidated, filterset, staging='ajcc', x_count="Filtered_Count", min_concentration=0.5, adjust_staging='early_late'):
    '''function plot_stage_vs_count() plots the stage vs the MMBIR count as a kde plot.'''
   
    #staging='ajcc'/'figo'/'tumor_grade', adjust='early_late'/'early_middle_late'

    df_consolidated = df_consolidated[df_consolidated["Sample_Type"].isin(filterset)]
    df_consolidated = df_consolidated[df_consolidated["Concentration"] >= min_concentration]

    print(df_consolidated.head())

    if staging == 'ajcc':
        col_name = "cases.0.diagnoses.0.ajcc_pathologic_stage"

    elif staging == 'figo':
        col_name = "cases.0.diagnoses.0.figo_stage"

    elif staging == 'tumor_grade':
        col_name = "cases.0.diagnoses.0.tumor_grade"

    else:
        print(f"Unknown staging input: {staging}. Please use 'ajcc' or 'figo'.")
        return
    
    df_consolidated.rename(columns={col_name: "Stage"}, inplace=True)
    
    #show the number of samples per stage
    print(df_consolidated.columns)
    print(df_consolidated["Stage"].value_counts())

    #if Stage X, remove it
    df_consolidated = df_consolidated[df_consolidated["Stage"] != "Stage X"]
    # if the stage is Stage I or Stage II, then set the stage to "Early"

    if adjust_staging=='early_late':

        print("Adjusting staging to early and late: Early = Stage IA, Stage IB, Stage IC, Stage I, Stage IIA, Stage IIB, Stage IIC; Late = everything else")
        df_consolidated["stage_adjusted"] = df_consolidated["Stage"].apply(lambda x: "Early" if x in ["Stage IA", "Stage IB", "Stage IC" "Stage I", "Stage IIA","Stage IIB", "Stage IIC"] else "Late")
    elif adjust_staging=='early_middle_late':

        print("Adjusting staging to early, middle and late: Early = Stage IA, Stage IB, Stage IC, Stage I; Middle = Stage IIA, Stage IIB, Stage IIC; Late = everything else")
        df_consolidated["stage_adjusted"] = df_consolidated["Stage"].apply(lambda x: "Early" if x in ["Stage IA", "Stage IB", "Stage IC", "Stage I"] else ("Middle" if x in ["Stage IIA","Stage IIB", "Stage IIC"] else "Late"))
    else:
        print(f"Unknown adjust_staging: {adjust_staging}")
        return
    
    #show the number of samples per stage
    print(df_consolidated["stage_adjusted"].value_counts())
    
    #show mean, median and standard deviation for each stage
    print("Mean, median and standard deviation for each stage")
    print(df_consolidated.groupby("stage_adjusted").mean()[x_count])
    print(df_consolidated.groupby("stage_adjusted").median()[x_count])
    print(df_consolidated.groupby("stage_adjusted").std()[x_count])

    #compare the mean of the early stage to the mean of the late stage statistically using a mann-whitney test
    print(stats.mannwhitneyu(df_consolidated[df_consolidated["stage_adjusted"]=="Early"][x_count], df_consolidated[df_consolidated["stage_adjusted"]=="Late"]["Raw_Count"]))
    print(stats.ttest_ind(df_consolidated[df_consolidated["stage_adjusted"]=="Early"][x_count], df_consolidated[df_consolidated["stage_adjusted"]=="Late"]["Raw_Count"]))

    sns.set_context("poster")
    sns.kdeplot(data=df_consolidated, x=x_count, hue="stage_adjusted", common_norm=False, common_grid=True)
    # sns.histplot(data=df_consolidated, x="Raw_Count", bins=200, hue="stage_adjusted", kde=True)
    plt.xlabel(f"{x_count} MMBIR Signature Count")
    plt.ylabel("Frequency")
    plt.show()

@fancy_status
def plot_stage_vs_concentration(df_consolidated, filterset, staging='ajcc', x_count="Filtered_Count", min_concentration=0.5, adjust_staging='early_late'):
    '''function plot_stage_vs_count() plots the stage vs concentration.
    The plot shows the concentration on the x-axis and the stage as a hist.
    you can adjust the staging to early, middle and late or early and late.
    staging='ajcc' or 'figo', adjust='early_late' or 'early_middle_late' '''

    df_consolidated = df_consolidated[df_consolidated["Sample_Type"].isin(filterset)]
    df_consolidated = df_consolidated[df_consolidated["Concentration"] >= min_concentration]

    if staging == 'ajcc':
        col_name="cases.0.diagnoses.0.ajcc_pathologic_stage"

    elif staging == 'figo':
        col_name="cases.0.diagnoses.0.figo_stage"

    elif staging == 'tumor_grade':
        col_name = "cases.0.diagnoses.0.tumor_grade"

    else:
        print(f"Unknown staging input: {staging}. Please use 'ajcc' or 'figo'.")
        return
    
    df_consolidated.rename(columns={col_name: "Stage"}, inplace=True)
    
    #show the number of samples per stage
    print(df_consolidated["Stage"].value_counts())

    #if Stage X, remove it
    df_consolidated = df_consolidated[df_consolidated["Stage"] != "Stage X"]

    if adjust_staging=='early_late':

        # if the stage is Stage I or Stage II, then set the stage to "Early"
        print("Adjusting staging to early and late: Early = Stage IA, Stage IB, Stage IC, Stage I, Stage IIA, Stage IIB, Stage IIC; Late = everything else")
        df_consolidated["stage_adjusted"] = df_consolidated["Stage"].apply(lambda x: "Early" if x in ["Stage IA", "Stage IB", "Stage IC" "Stage I", "Stage IIA","Stage IIB", "Stage IIC"] else "Late")

    elif adjust_staging=='early_middle_late':

        print("Adjusting staging to early, middle and late: Early = Stage IA, Stage IB, Stage IC, Stage I; Middle = Stage IIA, Stage IIB, Stage IIC; Late = everything else")
        df_consolidated["stage_adjusted"] = df_consolidated["Stage"].apply(lambda x: "Early" if x in ["Stage IA", "Stage IB", "Stage IC", "Stage I"] else ("Middle" if x in ["Stage IIA","Stage IIB", "Stage IIC"] else "Late"))
    
    else:
        print(f"Unknown adjust_staging: {adjust_staging}")
        return
    
    #show the number of samples per stage
    print(df_consolidated["stage_adjusted"].value_counts())
    
    #show mean, median and standard deviation for each stage
    print("Mean, median and standard deviation for each stage")
    print(df_consolidated.groupby("stage_adjusted").mean()[x_count])
    print(df_consolidated.groupby("stage_adjusted").median()[x_count])
    print(df_consolidated.groupby("stage_adjusted").std()[x_count])

    #show stage vs concentration
    sns.set_context("poster")
    sns.histplot(data=df_consolidated, x="Concentration", bins=200, hue="stage_adjusted", kde=True)
    plt.xlabel("Aliquot Concentration")
    plt.ylabel("Frequency")
    plt.show()

@fancy_status
def plot_ajcc_pathologic_n_vs_count(df_consolidated, filterset, x_count="Filtered_Count", min_concentration=0.5):
    '''function plot_ajcc_pathologic_n_vs_count() plots the ajcc pathologic n vs the MMBIR count.'''

    df_consolidated = df_consolidated[df_consolidated["Sample_Type"].isin(filterset)]
    df_consolidated = df_consolidated[df_consolidated["Concentration"] >= min_concentration]

    #rename columns
    df_consolidated.rename(columns={"cases.0.diagnoses.0.ajcc_pathologic_n": "StageN"}, inplace=True)

    #count how many samples of each stageN there are
    print(df_consolidated["StageN"].value_counts())

    #if Stage X, remove it
    df_consolidated = df_consolidated[df_consolidated["StageN"] != "NX"]

    #adjust the staging; group N0, N1, N1b, N1mi, N0 (i+) together as N0, N2+N3 together as N2+N3
    df_consolidated["stageN_adjusted"] = df_consolidated["StageN"].apply(lambda x: "N0" if x in ["N0", "N0 (i-)"] else ("N1" if x in ["N1","N1b","N1mi"] else ("N0 (i+)" if x in ["N0 (i+)"] else "N2+N3")))

    #show mean median and standard deviation for each stage
    print("Mean, median, standard deviation for each stage")
    print(df_consolidated.groupby("stageN_adjusted").mean()[x_count])
    print(df_consolidated.groupby("stageN_adjusted").median()[x_count])
    print(df_consolidated.groupby("stageN_adjusted").std()[x_count])

    #compare the mean of the early stage to the mean of the late stage statistically using a mann-whitney test
    print(stats.mannwhitneyu(df_consolidated[df_consolidated["stageN_adjusted"]=="N0"][x_count], df_consolidated[df_consolidated["stageN_adjusted"]=="N2+N3"]["Raw_Count"]))
    print(stats.ttest_ind(df_consolidated[df_consolidated["stageN_adjusted"]=="N0"][x_count], df_consolidated[df_consolidated["stageN_adjusted"]=="N2+N3"]["Raw_Count"]))

    sns.set_context("poster")
    sns.kdeplot(data=df_consolidated, x=x_count, hue="stageN_adjusted", common_norm=False, common_grid=True) #
    # sns.histplot(data=df_consolidated, x=x_count, bins=200, hue="stage_adjusted", kde=True)
    plt.xlabel("Raw MMBIR Signature Count")
    plt.ylabel("Frequency")
    plt.show()

@fancy_status
def plot_ajcc_pathologic_t_vs_count(df_consolidated, filterset, x_count="Filtered_Count", min_concentration=0.5):
    '''function plot_ajcc_pathologic_t_vs_count() plots the ajcc pathologic t vs the MMBIR count.'''

    df_consolidated = df_consolidated[df_consolidated["Sample_Type"].isin(filterset)]
    df_consolidated = df_consolidated[df_consolidated["Concentration"] >= min_concentration]

    df_consolidated.rename(columns={"cases.0.diagnoses.0.ajcc_pathologic_t": "StageT"}, inplace=True)

    #count how many samples of each stageN there are
    print("StageT counts:")
    print(df_consolidated["StageT"].value_counts())

    #df_consolidated = df_consolidated[df_consolidated["StageT"] != "T4d"]

    #bin the samples into T1, T2, T3&T4
    df_consolidated["stageT_adjusted"] = df_consolidated["StageT"].apply(lambda x: "T1" if x in ["T1", "T1c", "T1b"] else ("T2" if x in ["T2"] else ("T3&4"))) # if x in ["T3"] else "T4")))
    
    #count how many samples of each binned stageT there are
    print("StageT adjusted counts:")
    print(df_consolidated["stageT_adjusted"].value_counts())

    #show mean median and standard deviation for each stage
    print("Mean, median and standard deviation for each stage:")
    print(df_consolidated.groupby("stageT_adjusted").mean()[x_count])
    print(df_consolidated.groupby("stageT_adjusted").median()[x_count])
    print(df_consolidated.groupby("stageT_adjusted").std()[x_count])

    #compare the mean of the early stage to the mean of the late stage statistically using a mann-whitney test
    print(stats.mannwhitneyu(df_consolidated[df_consolidated["stageT_adjusted"]=="T1"][x_count], df_consolidated[df_consolidated["stageT_adjusted"]=="T3&4"]["Raw_Count"]))
    print(stats.ttest_ind(df_consolidated[df_consolidated["stageT_adjusted"]=="T1"][x_count], df_consolidated[df_consolidated["stageT_adjusted"]=="T3&4"]["Raw_Count"]))

    sns.set_context("poster")
    #sns.histplot(data=df_consolidated, x=x_count, bins=100, hue="stageT_adjusted", kde=True)
    sns.kdeplot(data=df_consolidated, x=x_count, hue="stageT_adjusted", common_norm=False, common_grid=True)
    plt.xlabel("Raw MMBIR Signature Count")
    plt.ylabel("Frequency")
    plt.show()

@fancy_status
def plot_age_vs_count_correlation(df_consolidated, count="Filtered_Count", min_concentration=0.5, method="spearman"):
    '''function plot_count_vs_age_correlation() plots the MMBIR count vs the age at collection.
    The plot shows the MMBIR count on the y-axis and the age at collection on the x-axis.'''

    df_consolidated = df_consolidated[df_consolidated["Concentration"] >= min_concentration]

    #if age at collection is not available, remove the sample
    df_consolidated = df_consolidated[df_consolidated["age_at_collection"].notna()]

    df_consolidated_tumor = df_consolidated[df_consolidated["Sample_Type"]=="Primary Tumor"]
    r_squared_tumor=df_consolidated_tumor.corr(method=method)["age_at_collection"][count]

    df_consolidated_blood = df_consolidated[df_consolidated["Sample_Type"]=="Blood Derived Normal"]
    r_squared_blood=df_consolidated_blood.corr(method=method)["age_at_collection"][count]

    print(f"R-squared for tumor: {r_squared_tumor}, R-squared for blood: {r_squared_blood}")

    #plot the correlation for both tumor and blood

    for df in [[df_consolidated_tumor, "Tumor"], [df_consolidated_blood, "Blood"]]:

        #calculate the linear regression between age_at_collection and count
        slope, intercept, r_value, p_value, std_err = stats.linregress(df[0]["age_at_collection"], df[0][count])
        print(f"Slope: {slope}, Intercept: {intercept}, R-squared: {r_value**2}, P-value: {p_value}")

        sns.set_context("poster")
        sns.scatterplot(x="age_at_collection", y=count, data=df[0], alpha=0.5)
        sns.regplot(x="age_at_collection",y=count, data=df[0], scatter=False, robust=True, color="orange")

        # add the pvalue and R-squared to the plot title, round the R-squared to 2 decimal places
        plt.title(f"R-squared: {r_value**2:.3f}, P-value: {p_value:.4f}, slope: {slope:.3f} for {df[1]}")

        # rename the x-axis and y-axis
        plt.xlabel("Age at sample collection (years)")
        plt.ylabel(f"MMBIR {count}")

        plt.show()

@fancy_status
def plot_age_vs_count_binned(df_consolidated, count="Filtered_Count", min_concentration=0.5):
    '''function plot_count_vs_age_binned() plots the MMBIR count vs the age at collection.
    The plot shows the MMBIR count on the y-axis and the age at collection on the x-axis.'''

    df_consolidated = df_consolidated[df_consolidated["Concentration"] >= min_concentration]

    #if age at collection is not available, remove the sample
    df_consolidated = df_consolidated[df_consolidated["age_at_collection"].notna()]

    #bin the samples into age bins
    # bins are 0-20, 20-40, 40-50, 50-60, 60-70, 70-80, 80+
    df_consolidated["age_bin"] = pd.cut(df_consolidated["age_at_collection"], bins=[0,20,40,50,60,70,80,200], labels=["0-20", "20-40", "40-50", "50-60", "60-70", "70-80", "80+"])

    #count how many samples of each age bin there are
    print("Age bin counts:")
    print(df_consolidated["age_bin"].value_counts())

    #show mean median and standard deviation for each age bin
    print("Mean, median and standard deviation for each age bin:")
    print(df_consolidated.groupby("age_bin").mean()[count])
    print(df_consolidated.groupby("age_bin").median()[count])
    print(df_consolidated.groupby("age_bin").std()[count])

    #compare the mean of the early age bin to the mean of the late age bin statistically using a mann-whitney test
    print(stats.mannwhitneyu(df_consolidated[df_consolidated["age_bin"]=="40-50"][count], df_consolidated[df_consolidated["age_bin"]=="80+"][count]))
    print(stats.ttest_ind(df_consolidated[df_consolidated["age_bin"]=="40-50"][count], df_consolidated[df_consolidated["age_bin"]=="80+"][count]))

    #plot the bins as a series of violin plots where the x-axis is the age bin and the y-axis is the count
    sns.set_context("poster")
    sns.violinplot(x="age_bin", y=count, data=df_consolidated, inner="quartile", scale="count", gridsize=1000)
    plt.xlabel("Age at sample collection (years)")
    plt.ylabel(f"MMBIR {count}")
    plt.show()

