!pip install datacache
!pip install pyfaidx
!pip install pyensembl
!pyensembl install --release 104 --species homo_sapiens

import pandas as pd
import ast
from pandas import DataFrame
import regex as re
from pyensembl import EnsemblRelease
from pyensembl.exon import Exon
import sys
from re import search

### PARAMS
census_dir='Cosmid.tsv'
mmb_df_input=sys.argv[1]

#Filtering

ref_complexity_filter=True
bir_complexity_filter=True
ref_homology_check_filter=True
bir_homology_check_filter=True


def unique(list1):
  # this function takes a list1 and returns a list2 with
  # unique elements of list 

  list_set = set(list1) 
  unique_list = (list(list_set))
  return unique_list

def wordsInSequence(sequence, treeLevel):
  # takes a genetic sequence and a treeLevel (word length) and returns
  # a list with words in that sequence

  length = len(sequence)
  max_possible_words_in_sequence = length-(treeLevel-1)
  wordList = []
  for i in range(max_possible_words_in_sequence):
    wordList.append(sequence[0+i:treeLevel+i])
  return wordList

def findComplexity(sequence, treeLevel, complexity_threshold):
  # takes a genetic sequence and calculates the linguistic complexity 
  # for that sequence for word length treeLevel

  wordList = wordsInSequence(sequence,treeLevel)
  wordList = unique(wordList)

  if len(sequence) > 4**treeLevel:
    complexity = len(wordList)/(4**treeLevel)
  else:
    complexity = len(wordList)/(len(sequence)-(treeLevel-1))

  if (complexity < complexity_threshold) & (complexity > 0.2):
    print("Complexity at tree level " + str(treeLevel) + " is " + str(complexity) + " for sequence: " + str(sequence))

  return ([complexity, treeLevel, sequence])

def movingWindow(sequence, treeLevel=8, chunkSize=20, complexity_threshold=0.2, showGraph = False):
  # takes a genetic sequence and calculates the linguistic complexity scores
  # for windows of size chunkSize for word lengths from 1 to treeLevel.
  # returns the lowest score window in format:  
  # [Boolean, [complexity, treeLevel, sequence]] where Boolean denotes if the
  # complexity score is lower (True) than complexity_threshold
  # optional: draw the complexity graph for the length of the sequence 
  # by setting showGraph = True

  complexityList=[]

  for i in range(1, treeLevel+1):
    for chunk in range(len(sequence)-(chunkSize-1)):

      chunkSeq = sequence[chunk:(chunk+chunkSize)]
      complexityList.append(findComplexity(chunkSeq, i, complexity_threshold))

  lowest=min(complexityList, key=lambda x: x[0])

  #print(complexityList)

  lowLvl=lowest[1]

  if showGraph:
    df = DataFrame(complexityList,columns=['complexity','tree','sequence'])
    df = df[['complexity','tree']]
    ax = df[df.tree == lowLvl].filter(items=["complexity"]).plot(fontsize=15, grid=True, figsize=(20,8))
    ax.axhline(complexity_threshold, color="red", linestyle="--")
    ax.axhline(1, color="green", linestyle="--")
    ax.set_ylim(ymin=0)
    ax.set_title(lowest[2])

  if lowest[0] < complexity_threshold:
    return([True, lowest[0], lowest[1]])
    #return([True, lowest])
  else:
    return([False, lowest[0], lowest[1]])
    #return([False, lowest])

def createLogFile(seqlist, filename, complexity_threshold = 0.2, chunkSize=20):
  # takes a list of sequences for analysis, the name of the output file,
  # complexity threshold, and windowsize (chunkSize), and outputs a log file
  # with lowest complexity score for each sequence
  outputf = open(filename, "w")
  truemmbir = 0
  for i in seqlist:
    seqScore = movingWindow(i, complexity_threshold = complexity_threshold, chunkSize = chunkSize)
    if seqScore[0] == False:
      truemmbir += 1
    outputf.write(str(seqScore[1])+"\n")
  print("sequences above threshold:", truemmbir)

  outputf.close()

##################
#Depreciated
def createDataFrameFromLogFile(logfile):
  # creates a pandas dataframe out of the complexity scores log file
  logfile = open(logfile, "r")
  lines = logfile.readlines()
  loglist=[]
  for line in lines:
    if "[" in line:
        loglist.append(line.strip("\n"))

  loglist2=[]
  for i in range(len(loglist)):
    loglist[i] = ast.literal_eval(loglist[i]) 
  print("total sequences:", len(loglist))

  df = DataFrame(loglist,columns=['complexity', 'level', 'sequence'])
  
  return df
################################################################################
# PART 2

def compl(base):
  if base == "A":
    return('T')
  elif base == "T":
    return('A')
  elif base == "G":
    return('C')
  elif base == "C":
    return('G')
  elif base == "-":
    return('-')

def rev_compl(seq):
  new_seq = ""
  for base in seq:
    new_base = compl(base)
    new_seq = new_base + new_seq
  return(new_seq)

def verifyImperfectHomology(ref, query, min_homology=0.8):
  mmbir = rev_compl(query)
  mmbir_errors = round(len(mmbir)*(1-min_homology))

  output_list = re.findall( '(' + mmbir + '){e<=' + str(mmbir_errors) + '}', ref)
  print(f'Found {len(output_list)} possible templates for {query}: {output_list}')

  if len(output_list) > 0:
    return(True)
  else:
    return(False)


def add_gene(row):

  data = EnsemblRelease(104)
  output = data.gene_names_at_locus(contig=row['chr'].strip("chr"), position=row["iBirStart"])
  #output = str(output)
  out_list = []
  for gname in output:
    if gname != "":
      out_list.append(gname)

  return(out_list) #output

def add_exon(row):

  data = EnsemblRelease(104)
  output = data.exons_at_locus(contig=row['chr'].strip("chr"), position=row["iBirStart"])
  #output = str(output)
  out_list = []
  for exon in output:
    if isinstance(exon, Exon):
      out_list.append(exon.exon_id)

  return(out_list) #output

def getMMBGenes(df_mmb_out,log=False):
  gene_list = df_mmb_out["genes"].tolist()

  outstr = ""
  gene_list_final=[]

  for inner_list in gene_list:
    for gene_name in inner_list:
      gene_list_final.append(gene_name)
      outstr+=f"{gene_name}\n"

  if log:
    print(outstr)
  gene_list_final = sorted(list(set(gene_list_final)))
  return(gene_list_final, outstr)


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


#Annotation of MMB Calls
df = pd.read_csv(mmb_df_input, sep="\t")

df[['ref_complexity_fail', 'ref_complexity_tree', 'ref_complexity_score']] = df.apply(lambda row: movingWindow(row.ref, complexity_threshold=0.2), axis=1, result_type ='expand')
df[['bir_complexity_fail', 'bir_complexity_tree', 'bir_complexity_score']] = df.apply(lambda row: movingWindow(row.bir, complexity_threshold=0.2), axis=1, result_type ='expand')
df['homology_check_ref'] = df.apply(lambda row: verifyImperfectHomology(row.ref, row.sBir), axis=1)
df['homology_check_bir'] = df.apply(lambda row: verifyImperfectHomology(row.bir, row.sBir), axis=1)
df["genes"] = df.apply(lambda row: add_gene(row), axis=1)
df["exones"] = df.apply(lambda row: add_exon(row), axis=1)


####
# filtering steps
filtered_df = df
if ref_complexity_filter:
  filtered_df = filtered_df[(filtered_df['ref_complexity_fail'] == False)]
if bir_complexity_filter:
  filtered_df = filtered_df[(filtered_df['bir_complexity_fail'] == False)]
if ref_homology_check_filter:
  filtered_df = filtered_df[(filtered_df['homology_check_ref'] == True)]
if bir_homology_check_filter:
  filtered_df = filtered_df[(filtered_df['homology_check_bir'] == True)]

filtered_df = filtered_df[filtered_df['exones'].map(lambda d: len(d)) > 0]


# finding cancer genes
gene_list_final, outstr = getMMBGenes(filtered_df,log=False)

print(f'the length of filtered gene list: {len(gene_list_final)}')
print(gene_list_final)
df_census = pd.read_csv(census_dir, sep='\t')

print("**************************************************\nCancer Genes:")
getCancerGeneNamesMMB(gene_list_final, df_census)
