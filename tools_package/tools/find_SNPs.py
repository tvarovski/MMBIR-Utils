import pandas as pd
import time
import sys

###REMINDER:
#we will create a multiprocessing pool to speed up the process (one process for each chromosome) in the future!
#we will also create a function to find the closest match SNP to each result

snps_df_path = "dbSNP_complex_indel8.csv"
results_df_path = sys.argv[1]

cancer_project = "OV"
filtered_status='raw' # 'raw' or 'filtered'


if results_df_path == "":
    print(f"no results_df_path given, using {cancer_project}\{filtered_status}_mmbir_results_master_df.csv")
    results_df_path = f"{cancer_project}\{filtered_status}_mmbir_results_master_df.csv"
else:
    print(f"job:{results_df_path}: started.")

def find_closest_snp(df_results, df_snps, jobid):
    #for each result in the df, find if there is a similar SNP in the df_snps 10bp around the result

    row_counter = 0
    time_start = time.time()
    time_now = time_start

    for index, row in df_results.iterrows():
        
        #print how many results have been processed
        if row_counter % 10000 == 1:
            print(f"job:{jobid}: processed {row_counter-1} results in {(time.time()-time_start):.4} seconds, {(time.time()-time_now):.4} seconds since last update.")
            time_now = time.time()
        row_counter += 1

        #get the chromosome and position of the result
        chrom = row['chr']
        pos = row['iBirStart']
        #get the SNPs that are within 10bp of the result
        df_snps_subset = df_snps[(df_snps['CHROM'] == chrom) & (df_snps['POS'] >= pos-10) & (df_snps['POS'] <= pos+10)].copy()

        #if there are any SNPs within 10bp of the result, pick the most similar SNP
        if len(df_snps_subset) > 0:
            #pick the SNP that is closest to the result
            df_snps_subset['distance'] = abs(df_snps_subset['POS'] - pos)
            df_snps_subset = df_snps_subset.sort_values(by='distance')
            #get the closest SNP
            closest_snp = df_snps_subset.iloc[0]
            #get the SNP ID
            snp_id = closest_snp['ID']
            #add the SNP ID to the df
            df_results.loc[index, 'SNP'] = snp_id
        else:
            df_results.loc[index, 'SNP'] = ''

    print(f"job:{jobid}: processed {row_counter} results in {time.time() - time_start} seconds.")
    return df_results

df = find_closest_snp(pd.read_csv(results_df_path), pd.read_csv(snps_df_path), results_df_path)

#save the df to a csv file (name based on the input file)
df.to_csv(f"{results_df_path}", index=False)
print(f"finished, saved to {results_df_path}")
