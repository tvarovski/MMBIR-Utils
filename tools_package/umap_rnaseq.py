import os
import logging
import umap
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import cancer_config as cfg
from tools import time_elapsed, fancy_status

#parameters#################

cancer = cfg.settings["TCGA-PROJECT"]
pickle_path = f"expression_data_{cancer}.pickle"
df_consolidated_path = f"consolidated_results_{cancer}.tsv"
username = cfg.settings["username"]

metadata_file = f"TCGA-{cancer}-WXS-BAM-metadata.tsv"
metadata_location = f"/Users/{username}/MMBIR_Databases/TCGA/{metadata_file}"

df_diff_expr_path = f"outputs/ttest_results_{cancer}_minconc0_bh_corrected.tsv"

sample_type = "Primary Tumor" #where to get MMBIR counts from
n_neighbors = 50
min_dist = 0.1
n_components = 2
metric = "canberra"

filter_transcriptome = True
min_pval = 0.05
fold_change_lower = 0.5
fold_change_upper = 2

scale_upper = 700
scale_lower = 50

show_plot = True

#############################

def filterByDiffExpr(df, df_diff_expr, min_pval=0.1, fold_change_lower=0.75, fold_change_upper=1.5):

    logging.basicConfig(level=logging.INFO)

    # take only the genes where the p-value is less than min_pval
    df_diff_expr = df_diff_expr[df_diff_expr["p-value"] < min_pval]
    logging.info(f"There are now {df_diff_expr.shape[0]} differentially expressed transcripts with p-value < {min_pval}.")

    # take only the genes where the fold change is outside fold_change_lower and fold_change_upper
    df_diff_expr = df_diff_expr[(df_diff_expr["fold-change"] < fold_change_lower) | (df_diff_expr["fold-change"] > fold_change_upper)]
    logging.info(f"There are now {df_diff_expr.shape[0]} differentially expressed transcripts with fold-change < {fold_change_lower} or > {fold_change_upper}.")

    #get the gene_id of the differentially expressed genes
    diff_expr_genes = df_diff_expr["gene_id"].tolist()

    #add columns to keep to the list
    diff_expr_genes.extend(['sample_name', 'sample_name_file', 'case_id'])

    #drop the columns from the df that are not in the list
    df = df[diff_expr_genes]

    return df

def addMMBToDF(df, df_consolidated, sample_type):

    logging.basicConfig(level=logging.INFO)

    # filter the df_consolidated dataframe to only include primary tumors
    df_consolidated = df_consolidated[df_consolidated["Sample_Type"] == sample_type]

    #group df_consolidated by Case_ID and get the highest value of Raw_Count and Filtered_Count
    grouping_dict = {"Raw_Count": "max", "Filtered_Count": "max"}
    df_consolidated = df_consolidated.groupby("Case_ID").agg(grouping_dict).reset_index()

    #add a column with Raw_Count to the df by merging on the Case_ID column (in df_consolidated) and the case_id column (in df)
    df = df.merge(df_consolidated[["Case_ID", "Raw_Count", "Filtered_Count"]], left_on="case_id", right_on="Case_ID")
    logging.info(f"Added the Raw_Count and Filtered_Count columns to the dataframe.")

    #remove the Case_ID column
    df = df.drop(columns=["Case_ID"])

    return df

def addMetadataToDF(df, df_metadata, columns):
    
    logging.basicConfig(level=logging.INFO)

    df = df.merge(df_metadata[["cases.0.case_id"]+columns], left_on="case_id", right_on="cases.0.case_id")
    logging.info(f"Added the {columns} columns to the dataframe.")

    #remove the "cases.0.case_id" column
    df = df.drop(columns=["cases.0.case_id"])
    
    return df

@fancy_status
@time_elapsed
def perform_UMAP(rnaseq_data, n_neighbors, min_dist, n_components, metric):

    logging.basicConfig(level=logging.INFO)

    logging.info("rna_seq_data shape:")
    logging.info(rnaseq_data.shape)

    # first, create a umap object
    reducer = umap.UMAP(        
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            n_components=n_components,
            metric=metric)

    logging.info("Scaling the data")
    scaled_rnaseq_data = StandardScaler().fit_transform(rnaseq_data)
    logging.info("Done scaling the data")

    # then, fit the umap object to the data
    logging.info("Fitting the umap object to the data")
    embedding = reducer.fit_transform(scaled_rnaseq_data)
    logging.info("Done fitting the umap object to the data")
    logging.info(f"The shape of the embedding is {embedding.shape}")
    logging.info("Done performing UMAP, returning the embedding")

    return reducer, embedding

def plot_UMAP_numerical(embedding, df, numerical_col, scale_upper=None, scale_lower=None, save_path=None, show_plot=True):

    logging.basicConfig(level=logging.INFO)

    # plot the umap projection, coloring by the Raw_Count column, set the highest value of color to 1000
    plt.figure(figsize=(12,12))
    # assign marker by figo_stage
    plt.scatter(embedding[:, 0], embedding[:, 1], c=df[numerical_col], s=16, alpha=0.8, cmap='viridis', vmin=scale_lower, vmax=scale_upper)
    plt.title(f'UMAP of RNAseq by {numerical_col}', fontsize=24, fontweight='bold')
    plt.xlabel("UMAP 1", fontsize=18)
    plt.ylabel("UMAP 2", fontsize=18)
    plt.axis('square')
    plt.colorbar(label=f"{numerical_col}")
    #save the figure
    if save_path is not None:
        plt.savefig(f"{save_path}.png", dpi=600)
        logging.info(f"Saved the figure to {save_path}")
    
    if show_plot:
        plt.show()
    else:
        plt.close()

def plot_UMAP_categorical(embedding, df, categorical_col, save_path=None, show_plot=True):

    logging.basicConfig(level=logging.INFO)

    color_labels = df[categorical_col].unique() # get the unique values of the column
    rgb_values = sns.color_palette("tab10", len(color_labels)) # create a list of rgb values
    color_map = dict(zip(color_labels, rgb_values)) # create a dictionary of the color labels and the rgb values
    logging.info(f"The color labels are {color_labels}")

    plt.figure(figsize=(12,12))
    # assign color by column_var which is a string, so it will be categorical
    plt.scatter(embedding[:, 0], embedding[:, 1], c=df[categorical_col].map(color_map), s=16, alpha=0.5)
    plt.axis('square')

    plt.title(f'UMAP of RNAseq by {categorical_col}', fontsize=24)
    plt.xlabel("UMAP 1", fontsize=18)
    plt.ylabel("UMAP 2", fontsize=18)

    # reuse the color_map for the legend
    patches = [mpatches.Patch(color=color, label=label) for label, color in color_map.items()]
    plt.legend(handles=patches, bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()

    if save_path is not None:
        plt.savefig(f"{save_path}", dpi=600)
        logging.info(f"Saved the figure to {save_path}")
    
    if show_plot:
        plt.show()
    else:
        plt.close()


# load the data from a pickle file, the consolidated results file, and the differentially expressed genes file
df = pd.read_pickle(pickle_path)
df_consolidated = pd.read_csv(df_consolidated_path, sep="\t")
df_diff_expr = pd.read_csv(df_diff_expr_path, sep="\t")
df_metadata = pd.read_csv(metadata_location, sep="\t")

if filter_transcriptome:

    df = filterByDiffExpr(df, df_diff_expr, min_pval=min_pval, fold_change_lower=fold_change_lower, fold_change_upper=fold_change_upper)

df = addMMBToDF(df, df_consolidated, sample_type)

# add the specified metadata to the dataframe
column_var_list_numerical = [   
                                "cases.0.samples.0.days_to_collection",
                                "cases.0.diagnoses.0.age_at_diagnosis",
                                "cases.0.demographic.days_to_birth"

                                ]
column_var_list_categorical = ["cases.0.disease_type",
                               "cases.0.diagnoses.0.figo_stage",
                               "cases.0.samples.0.portions.0.analytes.0.aliquots.0.concentration",
                               "cases.0.diagnoses.0.prior_malignancy",
                               "cases.0.diagnoses.0.prior_treatment",
                               "cases.0.diagnoses.0.tumor_grade",
                               "cases.0.samples.0.preservation_method",
                               "cases.0.samples.0.sample_type",
                               "cases.0.samples.0.sample_type_id",
                               "cases.0.samples.0.submitter_id",
                               "cases.0.submitter_id"
                               
                               ]
#"cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"

column_var_list=column_var_list_numerical + column_var_list_categorical

try:
    df = addMetadataToDF(df, df_metadata, column_var_list)

except Exception as e:
    logging.error(f"Error adding metadata: {column_var_list} to the dataframe: {e}")
    

#prepend "Raw_Count" and "Filtered_Count" to column_var_list_numerical
column_var_list_numerical = ["Raw_Count", "Filtered_Count"] + column_var_list_numerical
trailing_columns = 5 + len(column_var_list) # add the number of metadata columns to trailing_columns

#print last trailing_columns + 3 columns
print(df.iloc[:, -trailing_columns-3:])
print(f"There are {len(df.columns)} columns in the dataframe and {trailing_columns} are trailing")

# create a list of the columns that we want to use for the umap projection
#it is all columns except the last few columns which are not part of the data
rnaseq_data = df.iloc[:, :-trailing_columns]

# perform the UMAP projection
reducer, embedding = perform_UMAP(rnaseq_data, n_neighbors, min_dist, n_components, metric)


# create a folder in the outputs folder to store the plots
if not os.path.exists("outputs/plots/umap"):
    os.makedirs("outputs/plots/umap")


for column_var in column_var_list_numerical:

    if (column_var == "Raw_Count") or (column_var == "Filtered_Count"):

        scale_upper_pass = scale_upper
        scale_lower_pass = scale_lower

    else:
        scale_upper_pass = None
        scale_lower_pass = None

    try:
        plot_UMAP_numerical(embedding, 
                            df, 
                            column_var, 
                            scale_upper=scale_upper_pass, 
                            scale_lower=scale_lower_pass, 
                            save_path=f"outputs/plots/umap/umap_{cancer}_by_{column_var}.png", 
                            show_plot=show_plot)
    except:
        logging.info(f"Could not plot numerical_col {column_var}")
        pass


for column_var in column_var_list_categorical:

    if column_var == "cases.0.diagnoses.0.figo_stage":

        #replace the column_var names, specifically remove A, B, C, 1, 2
        df[column_var] = df[column_var].str.replace("A", "")
        df[column_var] = df[column_var].str.replace("B", "")
        df[column_var] = df[column_var].str.replace("C", "")
        df[column_var] = df[column_var].str.replace("1", "")
        df[column_var] = df[column_var].str.replace("2", "")
    try:
        plot_UMAP_categorical(embedding,
                            df, 
                            column_var, 
                            save_path=f"outputs/plots/umap/umap_{cancer}_by_{column_var}.png", 
                            show_plot=show_plot)
    except:
        logging.info(f"Could not plot categorical_col {column_var}")
        pass