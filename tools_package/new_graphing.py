import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_violin_for_all_TCGA_projects(path: str, data_type: str):
    
    if data_type == "raw":
        # Load the raw MMBIR data
        df_consolidated = pd.read_csv(path)
    elif data_type == "filtered":
        # Load the filtered MMBIR data
        df_consolidated = pd.read_csv(path, sep="\t")
    else:
        raise ValueError("data_type must be either 'raw' or 'filtered'")

    #columns:
    # Index(['chr', 'iBirStart', 'Consensus/cluster_number', 'iDepth', 'sBir',
    #        'sBirReversed', 'ref', 'bir', 'TEM', 'Microhomology_Insertion',
    #        'Microhomology_template', 'readStart', 'var16', 'var17', 'var18',
    #        'ref_complexity_fail', 'ref_complexity_score', 'ref_complexity_tree',
    #        'bir_complexity_fail', 'bir_complexity_score', 'bir_complexity_tree',
    #        'homology_check_ref', 'homology_check_bir', 'genes', 'exones',
    #        'average_read_length', 'case_id', 'days_to_birth', 'days_to_death',
    #        'vital_status', 'age_at_diagnosis', 'ajcc_pathologic_m',
    #        'ajcc_pathologic_n', 'ajcc_pathologic_stage', 'ajcc_pathologic_t',
    #        'tumor_grade', 'disease_type', 'project_id', 'days_to_collection',
    #        'aliquots.0.concentration', 'concentration', 'sample_type',
    #        'tumor_descriptor', 'file_name', 'id', 'mean_coverage',
    #        'proportion_reads_mapped', 'proportion_targets_no_coverage',
    #        'total_reads', 'SNP', 'df_name', 'figo_stage', 'prior_malignancy',
    #        'prior_treatment', 'preservation_method', 'sample_type_id',
    #        'submitter_id', 'tissue_type'],
    #       dtype='object')

    #group by the case_id, project_id, and sample_type columns and get the count of the number of rows in each group
    df_grouped = df_consolidated.groupby(['case_id', 'project_id', 'sample_type']).size().reset_index(name='counts')

    #create a violin plot for the counts within "Primary Tumor" and "Blood Derived Normal" for each project_id
    df_plot = df_grouped[df_grouped['sample_type'].isin(['Primary Tumor', 'Blood Derived Normal'])]

    #make plot wide to accommodate all project_ids
    plt.figure(figsize=(20, 10))

    sns.violinplot(
        x='project_id', y='counts', hue='sample_type', data=df_plot, split=True, inner="quartile",
        palette={"Primary Tumor": "lightblue", "Blood Derived Normal": "lightcoral"}, linewidth=1.25)

    plt.ylim(0, None)
    plt.xticks(rotation=90)
    plt.title(f"MMBIR Counts of Primary Tumor and Blood Derived Normal by Project ID ({data_type})")
    plt.xlabel("Project ID")
    plt.ylabel("Number of Raw MMBIR signatures")
    plt.legend(title="Sample Type")
    plt.savefig(f"counts_violin_plot_{data_type}.png", dpi=300, bbox_inches='tight')
    plt.show()
        
    print(f"Violin plot saved as counts_violin_plot_{data_type}.png")


data_type = "raw"  # "raw" or "filtered"
path = f"all_cancer_mmbir_{data_type}.csv"

plot_violin_for_all_TCGA_projects(path, data_type)

data_type = "filtered"  # "raw" or "filtered"
path = f"all_cancer_mmbir_{data_type}.csv"

plot_violin_for_all_TCGA_projects(path, data_type)



