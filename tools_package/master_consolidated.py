import pandas as pd


df_raw = pd.read_csv(f'all_mmbir_calls_raw.tsv', sep='\t')
df_filtered = pd.read_csv(f'all_mmbir_calls_filtered.tsv', sep='\t')
alliqut_concentration = 0.0

output_file = 'all_mmbir_calls_consolidated.tsv'

def add_df_filter(df: pd.DataFrame) -> pd.DataFrame:
    # Filter the DataFrame based on the specified conditions
    # Remove rows where 'ref_complexity_fail' is True or 'bir_complexity_fail' is True
    # Remove rows where 'homology_check_ref' is False or 'homology_check_bir' is False
    # Keep rows where 'exones' is not empty


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

    df = df[
        (df['ref_complexity_fail'] == False) &
        (df['bir_complexity_fail'] == False) &
        (df['homology_check_ref'] == True) &
        (df['homology_check_bir'] == True) &
        #if exones are an empty list, then remove the row too (but we want to keep it as a list)
        (~df['exones'].apply(lambda x: isinstance(x, list) and len(x) == 0)) &
        (df['exones'].apply(lambda x: str(x) != "[]")) &

        #if SNP is present (the format is a string), remove the row if it is empty or contains only whitespace
        (df['SNP'].apply(lambda x: str(x).strip() != ""))
    ]

    return df

def add_aliqut_filter(df: pd.DataFrame, alliqut_concentration: float = 0.0) -> pd.DataFrame:
    # Filter the DataFrame based on the specified conditions
    # Remove rows where 'aliquots.0.concentration' is less than the specified value
    df = df[df['aliquots.0.concentration'] >= alliqut_concentration]
    return df

df_filtered = add_df_filter(df_filtered)

df_raw = add_aliqut_filter(df_raw, alliqut_concentration=alliqut_concentration)
df_filtered = add_aliqut_filter(df_filtered, alliqut_concentration=alliqut_concentration)


#now we need to group by the columns to get the following columns:
# Project, Case_ID, Sample_Name, Sample_Type, Raw_Count, Filtered_Count

df_raw_grouped = df_raw.groupby(['project_id', 'case_id', 'sample_type']).size().reset_index(name='Raw_Count')
df_filtered_grouped = df_filtered.groupby(['project_id', 'case_id', 'sample_type']).size().reset_index(name='Filtered_Count')

df_merged = pd.merge(df_raw_grouped, df_filtered_grouped, on=['project_id', 'case_id', 'sample_type'], how='outer')

print("Null values in the merged dataframe:")
print(df_merged.isnull().sum())

#if there are nulls in Filtered_Count, then we need to fill them with 0
df_merged['Filtered_Count'] = df_merged['Filtered_Count'].fillna(0)
print("Filled null values in Filtered_Count with 0")
df_merged['Filtered_Count'] = df_merged['Filtered_Count'].astype(int)

# check if there are still null values in the merged dataframe
if df_merged.isnull().sum().sum() > 0:
    print("Warning: There are still null values in the merged dataframe")
else:
    print("No remaining null values in the dataframe")


#rename the columns to match the required format
df_merged.rename(columns={
    'project_id': 'Project', 
    'case_id': 'Case_ID', 
    'sample_type': 'Sample_Type'}, inplace=True)

# save the merged dataframe to a tsv file
df_merged.to_csv(output_file, sep='\t', index=False)
print(f"Merged dataframe saved to {output_file}")
