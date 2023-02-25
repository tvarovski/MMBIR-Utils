import os
import pandas as pd
import logging
from tools import getCasesAboveMMBThreshold, fancy_status, time_elapsed, findThresholdCases

logger = logging.getLogger(__name__)
logging.basicConfig()

@fancy_status
@time_elapsed
def masked_snv_mv(df_metadata, snvs_loc=""):

  output = os.system(f"ls -d1 *-*-*-*-* > cases.txt")
  print(f"cases exited with a code {output}")
  cases_file = open("cases.txt", 'r')
  cases = cases_file.readlines()
  cases_file.close()


  for caseID in cases:
      
    caseID=caseID.strip()
    case_snvs = df_metadata[df_metadata["cases.0.case_id"] == caseID]
    case_files = case_snvs[["file_name", "id"]].values.tolist()
    for file_name, id in case_files:
      output = os.system(f"cp {snvs_loc}/{id}/{file_name} {caseID}")
      if output == 1:
        print(f"there is a problem with {snvs_loc}/{id}/{file_name} {caseID}")
      output = os.system(f"yes n | gunzip {caseID}/{file_name}")

def extractAffectedGeneNames(maf_path):

    try:
        df_masked_snvs=pd.read_csv(maf_path, sep="\t", comment='#')
    except OSError:
        logging.info(f"Couldnt open the MAF file: {maf_path} Is it there?")
        return None
    try:
        # remove silent mutations
        df_masked_snvs = df_masked_snvs[df_masked_snvs["Variant_Classification"] != "Silent"]
        return_gene_list = df_masked_snvs["Hugo_Symbol"].values.tolist()
        return return_gene_list
    except ValueError:
        logging.info(f"Couldn't parse the MAF file: {maf_path}")
        return None

def findCurrentCases():
    output = os.system(f"ls -d1 *-*-*-*-* > cases.txt")
    print(f"cases exited with a code {output}")
    cases_file = open("cases.txt", 'r')
    cases = cases_file.readlines()
    cases_file.close()
    return cases

def getSNV_frequencies(cases, df_metadata_maf, threshold_mmb_cases_low_df, threshold_mmb_cases_high_df):
    import numpy

    #initiate variables for the high and low mmbir cohorts
    high_mmbir_cases=0
    high_mmbir_snv_genes=[]
    low_mmbir_cases=0
    low_mmbir_snv_genes=[]

    for caseID in cases:

        caseID=caseID.strip()

        case_snvs = df_metadata_maf[df_metadata_maf["cases.0.case_id"] == caseID]
        case_files = case_snvs[["file_name", "id"]].values.tolist()

        #get all the snv affected genes (for now just pick the first maf file comming in)
        for file_name, id in case_files:
            maf_path = f"{caseID}/{file_name}"
            maf_path = maf_path.strip(".gz")
            maf_genes = extractAffectedGeneNames(maf_path)


        #check if case is high or low mmbir
        if caseID in threshold_mmb_cases_high_df["Case_ID"].values.tolist():

            #if yes -> HIGH MMBIR, count as high mmbir sample, add to high mmbir cohort
            high_mmbir_cases+=1
            high_mmbir_snv_genes.extend(maf_genes)

        elif caseID in threshold_mmb_cases_low_df["Case_ID"].values.tolist():

            #if no -> LOW MMBIR, count as low mmbir sample, add to low mmbir cohort
            low_mmbir_cases+=1
            low_mmbir_snv_genes.extend(maf_genes)

        else:
            #if neither -> set is ambiguous, do nothing
            pass
    
    #get the average number of affected genes per case
    avg_high_mmbir_snv_genes = len(high_mmbir_snv_genes)/high_mmbir_cases
    avg_low_mmbir_snv_genes = len(low_mmbir_snv_genes)/low_mmbir_cases

    print(f"average number of affected genes per case in HIGH: {avg_high_mmbir_snv_genes} and LOW: {avg_low_mmbir_snv_genes}")


    #calculate the frequency distribution of the genes in the high, low and all sets
    fdist_HIGH=dict(zip(*numpy.unique(high_mmbir_snv_genes, return_counts=True)))
    fdist_HIGH = {k: v / high_mmbir_cases for k, v in fdist_HIGH.items()}
    sorted_fdist_HIGH = sorted(fdist_HIGH.items(), key=lambda x: x[1], reverse=True)

    fdist_LOW=dict(zip(*numpy.unique(low_mmbir_snv_genes, return_counts=True)))
    fdist_LOW = {k: v / low_mmbir_cases for k, v in fdist_LOW.items()}
    sorted_fdist_LOW = sorted(fdist_LOW.items(), key=lambda x: x[1], reverse=True)

    fdist_all=dict(zip(*numpy.unique(high_mmbir_snv_genes+low_mmbir_snv_genes, return_counts=True)))
    fdist_all = {k: v / (high_mmbir_cases+low_mmbir_cases) for k, v in fdist_all.items()}
    sorted_fdist_all = sorted(fdist_all.items(), key=lambda x: x[1], reverse=True)


    print(f"found HIGH: {high_mmbir_cases} and LOW: {low_mmbir_cases} cases")
    print("HIGH")
    print("Hugo,high-set,low-set,all-set")

    #create an empty dataframe to store the results
    df = pd.DataFrame(columns=["Hugo", "high-set", "low-set", "all-set"])

    for i in sorted_fdist_HIGH:

        try:
            compl_low=fdist_LOW[i[0]]

        except:
            logging.info(f"No corresponding record found for {i[0]} in the LOW-MMBIR Set")
            compl_low=0

        compl_all=fdist_all[i[0]]

        print(f"{i[0]},{i[1]},{compl_low},{compl_all}")

        # add the results to the dataframe, concat
        df = pd.concat([df, pd.DataFrame([[i[0], i[1], compl_low, compl_all]], columns=["Hugo", "high-set", "low-set", "all-set"])])
        
    print("LOW")
    for i in sorted_fdist_LOW:

        try:
            compl_high=fdist_HIGH[i[0]]

        except:
            logging.info(f"No corresponding record found for {i[0]} in the HIGH-MMBIR Set")
            compl_high=0

        compl_all=fdist_all[i[0]]

        print(f"{i[0]},{compl_high},{i[1]},{compl_all}")

        # add the results to the dataframe, concat
        df = pd.concat([df, pd.DataFrame([[i[0], compl_high, i[1], compl_all]], columns=["Hugo", "high-set", "low-set", "all-set"])])

    # remove duplicate rows from the dataframe
    df = df.drop_duplicates()

    #sort the dataframe by the high-set column, largest to smallest
    df = df.sort_values(by="high-set", ascending=False)

    #check if the Hugo column contains unique values
    if df["Hugo"].is_unique:
        print("Hugo column contains unique values")
    else:
        print("Hugo column contains duplicate values")

    # save the dataframe to a csv file

    return df

@fancy_status
@time_elapsed
def performSNVanalysis(params):
    consolidated_results_path = params["consolidated_results_path"]
    df_metadata = params["df_metadata"]
    df_metadata_maf = params["df_metadata_maf"]
    MMBIR_THRESHOLD_LOW = params["MMBIR_THRESHOLD_LOW"]
    MMBIR_THRESHOLD_HIGH = params["MMBIR_THRESHOLD_HIGH"]
    fraction_high= params["mmbir_fraction_high"]
    fraction_low= params["mmbir_fraction_low"]
    outputs_path = params["outputs_path"]

    min_concentration = params["min_concentration"]

    #old method of getting high+low mmbir cases
    #threshold_mmb_cases_high_df = getCasesAboveMMBThreshold(consolidated_results_path, df_metadata, MMBIR_THRESHOLD_HIGH, min_concentration=min_concentration)
    #threshold_mmb_cases_low_df = getCasesAboveMMBThreshold(consolidated_results_path, df_metadata, MMBIR_THRESHOLD_LOW, below=True, min_concentration=min_concentration)

    #print(f"found {len(threshold_mmb_cases_high_df)} cases above threshold {MMBIR_THRESHOLD_HIGH}")
    #print(f"found {len(threshold_mmb_cases_low_df)} cases below threshold {MMBIR_THRESHOLD_LOW}")

    #new method of getting high+low mmbir cases
    df_consolidated = pd.read_csv(consolidated_results_path, sep="\t")

    #default behavior is to filter out cases based on Raw_Counts not Filtered_Counts, but this can be changed by setting filtered=True
    threshold_mmb_cases_high_df, threshold_mmb_cases_low_df = findThresholdCases(df_consolidated, df_metadata, fraction_high=fraction_high, fraction_low=fraction_low, min_concentration=min_concentration, filtered=False)

    cases = findCurrentCases()

    df = getSNV_frequencies(cases, df_metadata_maf, threshold_mmb_cases_low_df, threshold_mmb_cases_high_df)
    df.to_csv(f"{outputs_path}/mmbir_snv_fraction_low_{fraction_low}_high_{fraction_high}.csv", index=False)