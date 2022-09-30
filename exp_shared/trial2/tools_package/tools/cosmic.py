from re import search
import sys
import pandas as pd
import ast

### PARAMS

ref_complexity_filter=True
bir_complexity_filter=True
ref_homology_check_filter=True
bir_homology_check_filter=True
exones_only=True

#ref_complexity_filter=False
#bir_complexity_filter=False
#ref_homology_check_filter=False
#bir_homology_check_filter=False
#exones_only=False


census_dir='tools/Cosmid.tsv'
mmb_df_input=sys.argv[1]

### END PARAMS



def getCancerGeneNamesMMB(genes_list, df_census):

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



def getMMBGenes(df_mmb_out):
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


def df_filter(df):

  filtered_df = df
  if ref_complexity_filter:
    filtered_df = filtered_df[(filtered_df['ref_complexity_fail'] == False)]
  if bir_complexity_filter:
    filtered_df = filtered_df[(filtered_df['bir_complexity_fail'] == False)]
  if ref_homology_check_filter:
    filtered_df = filtered_df[(filtered_df['homology_check_ref'] == True)]
  if bir_homology_check_filter:
    filtered_df = filtered_df[(filtered_df['homology_check_bir'] == True)]
  if exones_only:
    filtered_df = filtered_df[filtered_df['exones'].map(lambda d: len(d)) > 0]
  
  return filtered_df

#Filter of MMB Calls
df = pd.read_csv(mmb_df_input, sep="\t")
filtered_df=df_filter(df)
gene_list_final, outstr = getMMBGenes(filtered_df)


# finding cancer genes

print(f'the length of filtered gene list: {len(gene_list_final)}')
print(gene_list_final)
#print(outstr)

print("**************************************************\nCancer Genes:")

df_census = pd.read_csv(census_dir, sep='\t')
getCancerGeneNamesMMB(gene_list_final, df_census)
