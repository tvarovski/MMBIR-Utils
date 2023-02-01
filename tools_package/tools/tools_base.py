import os
import pandas as pd

def groupCases(df_metadata):
  import shutil 
  mypath = os.path.abspath(os.getcwd())

  onlydirs = [d for d in os.listdir(mypath) if os.path.isdir(os.path.join(mypath, d))]

  print(f"Found {len(onlydirs)} directories: {onlydirs}")

  for current_name in onlydirs:

    current_file_entry = df_metadata[df_metadata["file_name"] == f"{current_name}.bam"]
    current_case_id = current_file_entry["cases.0.case_id"].tolist()
    if len(current_case_id) == 0:
      continue
    if len(current_case_id) != 1:
      print("Retreived unexpected number of case_id for sample metadata. Exitting...")
      print(f"Case IDs: {current_case_id}")
      exit()
    current_case_id_samples = df_metadata[df_metadata["cases.0.case_id"] == current_case_id[0]]

    current_case_id_samples_list = current_case_id_samples["file_name"].tolist()

    for sample in current_case_id_samples_list:
      current_case_id_path = os.path.join(mypath, current_case_id[0])
      if not os.path.exists(current_case_id_path):
        os.makedirs(current_case_id_path)

      if os.path.exists(os.path.join(mypath, sample[:-4])):
        try:
          shutil.move(sample[:-4], current_case_id_path)
          print(f"Moved {sample[:-4]} to {current_case_id_path}")
        except:
          print("directory exist, couldn't move it")

def parseOutputs(df_metadata, consolidated_results_name="consolidated_results.tsv"):

  df_metadata = df_metadata.sort_values(by=["cases.0.case_id"])

  #create a table for raw and filtered numbers of MMBIR events (headers and numbers)
  output = os.system(f"echo 'sample	mmbir_events' > mmbir_counts_raw.tsv")
  output = os.system(f"./count_events.sh outputs/raw >> mmbir_counts_raw.tsv")

  output = os.system(f"echo 'sample	mmbir_events' > mmbir_counts_filtered.tsv")
  output = os.system(f"./count_events.sh outputs/filtered >> mmbir_counts_filtered.tsv")

  raw_df = pd.read_csv("mmbir_counts_raw.tsv", sep="\t")
  raw_df["sample"]=raw_df["sample"].str[:-8]
  filtered_df = pd.read_csv("mmbir_counts_filtered.tsv", sep="\t")
  filtered_df["sample"]=filtered_df["sample"].str[:-29]

  output = os.system(f"ls -d1 *-*-*-*-* > cases.txt")
  print(f"cases exited with a code {output}")
  cases_file = open("cases.txt", 'r')
  cases = cases_file.readlines()
  cases_file.close()


  column_names = ["Case_ID", "Sample_Name", "Sample_Type", "Raw_Count", "Filtered_Count"]
  output_df=pd.DataFrame(columns = column_names)

  for caseID in cases:
    caseID=caseID.strip()
    case_files = df_metadata[df_metadata["cases.0.case_id"] == caseID]
    case_files = case_files[["file_name", "cases.0.samples.0.sample_type"]]
    output = os.system(f"ls --hide=*.* -1 {caseID} > {caseID}/ind_files.txt")

    samples_file = open(f"{caseID}/ind_files.txt", 'r')
    samples = samples_file.readlines()
    samples_file.close()
    
    parsing_list=[]
    for sample in samples:
      sample=f"{sample.strip()}.bam"
      df_sample = case_files[case_files["file_name"] == sample]
      sample_list = df_sample[["file_name", "cases.0.samples.0.sample_type"]].values.tolist()
      try:
        parsing_list.append(sample_list[0])
      except IndexError:
        print(f"Couldn't read data in sample list: {sample_list} of CaseID: {caseID}")
        exit()

    for i in parsing_list:

      raw_count_list = raw_df[raw_df["sample"]==i[0][:-4]]
      raw_count_list = raw_count_list["mmbir_events"].tolist()
      filtered_count_list = filtered_df[filtered_df["sample"]==i[0][:-4]]
      filtered_count_list = filtered_count_list["mmbir_events"].tolist()
      if len(filtered_count_list) == 0:
        filtered_count_list = 0
      append_me=pd.DataFrame({"Case_ID": caseID,
                              "Sample_Name": i[0],
                              "Sample_Type": i[1],
                              "Raw_Count": raw_count_list,
                              "Filtered_Count": filtered_count_list})

      output_df = pd.concat([append_me,output_df.loc[:]]).reset_index(drop=True)
    output = os.system(f"rm {caseID}/ind_files.txt")
  output = os.system(f"rm cases.txt")

  print(output_df)
  output_df.to_csv(consolidated_results_name, sep="\t", index=False)
  
def masked_snv_mv(df_metadata):

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

def getCasesAboveMMBThreshold(consolidated_results_path, min_MMBIR_events, below=False, min_concentration=0):

    df_consolidated=pd.read_csv(consolidated_results_path, sep="\t")

    df_consolidated=pd.merge(df_consolidated, df_sample_metadata, left_on="Sample_Name", right_on="file_name")

    df_consolidated["age_at_collection"] = df_consolidated["cases.0.diagnoses.0.age_at_diagnosis"] + df_consolidated["cases.0.samples.0.days_to_collection"]
    df_consolidated.rename(columns={"cases.0.samples.0.portions.0.analytes.0.aliquots.0.concentration": "Concentration"}, inplace=True)
    df_consolidated=df_consolidated[df_consolidated["Concentration"] >= min_concentration]


    agg_dict={"Raw_Count": ['min', 'max'],
              "Filtered_Count": ['min', 'max']}
    df_agg = df_consolidated.groupby("Case_ID").agg(agg_dict).reset_index()

    df_agg.columns = ['_'.join(col).strip() for col in df_agg.columns.values]

    if below:
        df_agg = df_agg[df_agg["Filtered_Count_max"] < min_MMBIR_events]

    else:
        df_agg = df_agg[df_agg["Filtered_Count_max"] >= min_MMBIR_events]


    logging.info(f"The length is: {len(df_agg)}")

    df_agg = df_agg.sort_values("Filtered_Count_max",ascending=False)

    return df_agg

def returnGenes(file, high_mmbir, exones_only=False, homology_ref=False, homology_bir=False, complexity_ref=False, complexity_bir=False):

    df = pd.read_csv(file, sep='\t')

    if exones_only:
        df = df[df["exones"].str.len() != 2]
        #print(df)
    if homology_ref:
        df = df[df["homology_check_ref"] == True]
    if homology_bir:
        df = df[df["homology_check_bir"] == True]
    if complexity_ref:
        df = df[df["ref_complexity_fail"] == False]
    if complexity_bir:
        df = df[df["bir_complexity_fail"] == False]

    # read genes column of the df and return list of genes
    genes = df["genes"].tolist()
    # get unique genes
    genes = list(set(genes))

    #if mmbir is high, return each gene with a high flag
    if high_mmbir == True:
        for index, gene in enumerate(genes):
            genes[index] = (gene, True)
    
    #if mmbir is low, return each gene with a low flag
    elif high_mmbir == False:
        for index, gene in enumerate(genes):
            genes[index] = (gene, False)
    
    else:
        raise ValueError("high_mmbir must be True or False")

    return genes

def countGenes(genes):
    gene_count_high = {}
    gene_count_low = {}
    for gene, high_mmbir in genes:
        if high_mmbir == True:
            if gene in gene_count_high:
                gene_count_high[gene] += 1
            else:
                gene_count_high[gene] = 1

        if high_mmbir == False:
            if gene in gene_count_low:
                gene_count_low[gene] += 1
            else:
                gene_count_low[gene] = 1
    return gene_count_high, gene_count_low

def get_case_bam_files(case_id, df_metadata):
    #this function returns a list of bam files for a given case_id
    df_case = df_metadata[df_metadata["cases.0.case_id"]==case_id]
    bam_files = df_case["file_name"].tolist()
    return bam_files

def get_bam_file_metadata(bam_file, df_metadata):
    #this function returns metadata associated with a given bam file
    df_bam = df_metadata[df_metadata["file_name"]==bam_file]
    #return a dictionary of metadata
    return df_bam

def check_for_missing_bams(df_metadata):
    #go through all the cases and create a list of all the files
    cases_list = []
    for index, row in df_metadata.iterrows():
        case_id = row['cases.0.case_id']
        cases_list.append(case_id)

    #remove duplicates
    cases_list = list(dict.fromkeys(cases_list))

    #check if the file exists as a directory
    found_cases=0
    found_bams=0
    missing_bams = []
    missing_cases = []
    for case_id in cases_list:
        #get current directory
        current_dir = os.getcwd()
        case_dir = f"{current_dir}/{case_id}"
        if not os.path.exists(case_dir):
            print(f"WARNING: The directory: {case_dir} does not exist")
            missing_cases.append(case_dir)
        else:
            found_cases+=1
            #enter the directory
            os.chdir(case_dir)
            #check what directories are in the directory
            dir_list = os.listdir()
            
            #get the bam files for the case
            case_bam_files = get_case_bam_files(case_id, df_metadata)
            #check if the bam files are in the directory
            for bam_file in case_bam_files:
                if bam_file.strip(".bam") not in dir_list:
                    print(f"WARNING: The bam file: {bam_file} does not exist in: {case_dir}")
                    missing_bams.append(bam_file)
                else:
                    found_bams+=1
            #go back to the current directory
            os.chdir(current_dir)
    print(f"Found {found_cases} cases out of {len(cases_list)}")
    print(f"Found {found_bams} bam files out of {len(df_metadata)}")
    print("\nMissing bams:\n")
    print(missing_bams)
    print("Missing cases:\n")
    print(missing_cases)
    return([missing_bams, missing_cases])

def create_missing_bams_manifest(missing_files, manifest_location, missing_manifest_output_name):

    if len(missing_files[0]) == 0 & len(missing_files[0]) == 0:
        print("All good! All files present!")
        return
    else:
        # create an empty, placeholder dataframe
        df_out = pd.DataFrame()
        #open the manifest file
        df = pd.read_csv(manifest_location, sep="\t")
        if len(missing_files[0]) > 0:
            df_missing_bams = df[df["filename"].isin(missing_files[0])]
            df_out = pd.concat([df_out, df_missing_bams], axis=0)
            
        if len(missing_files[1]) > 0:
            #open the metadata file to retreive filenames associated with missing caseIDs
            df_metadata = pd.read_csv(metadata_location, sep="\t")
            df_metadata_missing_cases = df_metadata[df_metadata["cases.0.case_id"].isin(missing_files[1])]

            #get the manifest rows containing only files from missing cases
            df_missing_cases = df[df["filename"].isin(df_metadata_missing_cases["file_name"])]
            df_out = pd.concat([df_out, df_missing_cases], axis=0)
        
        #save the df_out to the new output_file
        df_out.to_csv(missing_manifest_output_name, sep="\t", index=False)
        return(df_out)

def create_mmbir_results_master_df(df_metadata, filtered=False, log=False):
    # create a master dataframe for raw mmbir results

    results_master_df = pd.DataFrame()
    if filtered:
        result_dir="filtered"
    else:
        result_dir="raw"
    #get current directory
    current_dir = os.getcwd()
    os.chdir(f"{current_dir}/outputs/{result_dir}")
    #list the files in the directory
    dir_list = os.listdir()
    #go through the files in the directory
    for file_name in dir_list:

        file_name = file_name.split("/")[-1]
        bam_name = file_name.replace("_raw.txt", "")
        bam_name = bam_name.replace("_allchrmmbir_filtered_all.txt", "")
        bam_name += ".bam"

        bam_metadata = get_bam_file_metadata(bam_name, df_metadata)
        #print(bam_metadata)

        if len(bam_metadata) == 0:
            print(f"WARNING: No metadata found for: {bam_name}, skipping")
            continue
        elif len(bam_metadata) > 1:
            print(f"WARNING: Multiple metadata found for: {bam_name}, skipping")
            continue

        #get column names of the bam_metadata dataframe
        column_names = bam_metadata.columns.tolist()

        #read in the filtered mmbir results
        filtered_mmbir_results_df = pd.read_csv(file_name, sep="\t")

        #add the bam_metadata to the filtered_mmbir_results_df, rename the column names for brevity
        for column_name in column_names:
            new_column_name = column_name.replace("cases.0.", "")
            new_column_name = new_column_name.replace("demographic.", "")
            new_column_name = new_column_name.replace("diagnoses.0.", "")
            new_column_name = new_column_name.replace("project.", "")
            new_column_name = new_column_name.replace("samples.0.", "")
            new_column_name = new_column_name.replace("portions.0.", "")
            new_column_name = new_column_name.replace("analytes.0.", "")

            #add the bam metadata to the dataframe
            filtered_mmbir_results_df[new_column_name] = bam_metadata[column_name].tolist()[0]
        
        #concat the filtered mmbir results to the master dataframe
        results_master_df = pd.concat([results_master_df, filtered_mmbir_results_df], axis=0)
        if log:
            print(f"Finished processing: {bam_name}")
    
    #go back to the current directory
    os.chdir(current_dir)
    return results_master_df

def extractAffectedGeneNames(maf_path):
    import logging

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

def createFullCancerTable(params):
    #if not working, uncomment the os library lines

    #unwrap the dictionary of params
    df_metadata = params["df_metadata"]
    manifest_location = params["manifest_location"]
    missing_manifest_output_name = params["missing_manifest_output_name"]
    log = params["log"]
    output_raw_name = params["output_raw_name"]
    output_filtered_name = params["output_filtered_name"]

    #my_dir=os.getcwd()
    #check for missing bams
    missing_files = check_for_missing_bams(df_metadata)
    #pass the list, if missing create manifest for redownload
    create_missing_bams_manifest(missing_files, manifest_location, missing_manifest_output_name)
    #create the master dataframe for the raw mmbir results
    raw_mmbir_results_master_df = create_mmbir_results_master_df(df_metadata, filtered=False, log=log)
    print(f"Finished creating the master dataframe for the raw mmbir results")
    #create the master dataframe for the filtered mmbir results
    filtered_mmbir_results_master_df = create_mmbir_results_master_df(df_metadata, filtered=True, log=log)
    print(f"Finished creating the master dataframe for the filtered mmbir results")

    #os.chdir(my_dir)
    #save the master dataframes to csv
    raw_mmbir_results_master_df.to_csv(output_raw_name, index=False)
    filtered_mmbir_results_master_df.to_csv(output_filtered_name, index=False)
    print(f"Finished saving the master dataframes to csv")

#Used by findCosmicGenes
def getCancerGeneNamesMMB(genes_list, df_census):

    from re import search

    df_census["Synonyms"] = df_census["Synonyms"].astype(str)

    for i in genes_list:

        for ind in df_census.index:

            if search("[\\^, ]"+i+"[, $]", df_census['Gene Symbol'][ind]):
                
                print(i,df_census['Gene Symbol'][ind], df_census['Synonyms'][ind],
                      df_census['Name'][ind], df_census['Tier'][ind])
            
            else:

                if search("[\\^, ]"+i+"[, $]", df_census['Synonyms'][ind]):

                    print(i,df_census['Gene Symbol'][ind], df_census['Synonyms'][ind],
                          df_census['Name'][ind], df_census['Tier'][ind])
                else:

                    if df_census['Gene Symbol'][ind] == i:

                        print(i,df_census['Gene Symbol'][ind], df_census['Synonyms'][ind],
                            df_census['Name'][ind], df_census['Tier'][ind])

#Used by findCosmicGenes
def getMMBGenes(df_mmb_out):
    import ast

    gene_list = df_mmb_out["genes"].tolist()
    outstr = ""
    gene_list_final=[]

    for gene_name in gene_list:

        gene_names = ast.literal_eval(gene_name)

        for gene_name in gene_names:

            gene_list_final.append(gene_name)
            outstr+=f"{gene_name}\n"

    gene_list_final = sorted(list(set(gene_list_final)))

    return(gene_list_final, outstr)

#Used by findCosmicGenes
def df_filter(df, filter_dict):

    filtered_df = df

    if filter_dict["ref_complexity_filter"]:
        filtered_df = filtered_df[(filtered_df['ref_complexity_fail'] == False)]

    if filter_dict["bir_complexity_filter"]:
        filtered_df = filtered_df[(filtered_df['bir_complexity_fail'] == False)]

    if filter_dict["ref_homology_check_filter"]:
        filtered_df = filtered_df[(filtered_df['homology_check_ref'] == True)]

    if filter_dict["bir_homology_check_filter"]:
        filtered_df = filtered_df[(filtered_df['homology_check_bir'] == True)]

    if filter_dict["exones_only"]:
        filtered_df = filtered_df[filtered_df['exones'].map(lambda d: len(d)) > 0]
    
    return filtered_df

#Used by findCosmicGenes
def findCosmicGenes(mmb_df_input_path, filter_dict):

    #Filter of MMB Calls
    df = pd.read_csv(mmb_df_input_path, sep="\t")
    filtered_df=df_filter(df, filter_dict)
    gene_list_final, outstr = getMMBGenes(filtered_df)

    # finding cancer genes
    print(f"the length of filtered MMBIR list: {len(filtered_df)}")
    print(f'the length of filtered gene list: {len(gene_list_final)}')
    print(gene_list_final)
    print("**************************************************\nCancer Genes:")

    df_census = pd.read_csv(census_dir, sep='\t')
    getCancerGeneNamesMMB(gene_list_final, df_census)
