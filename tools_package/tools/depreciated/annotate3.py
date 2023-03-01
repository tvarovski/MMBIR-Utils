# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

import pandas as pd
import ast
from pandas import DataFrame
from pyensembl import EnsemblRelease
from pyensembl.exon import Exon
import sys
import numpy as np
import regex as re

#PARAMS:

filter = sys.argv[3]
if filter == "False":
  f_name = "final_bir_locs.txt"
elif filter == "True":
  f_name = "non_commons_fromAll.txt"
else:
  print("argument not recognized. Exitting")
  exit()


path=sys.argv[1]
#End PARAMS

def createChrList(chrNum):
  chrList = []
  for i in range(chrNum):
    chrList.append(f"chr{i+1}")
  return(chrList)

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
  if mmbir_errors > 10:
    mmbir_errors = 10
  output_list = re.findall( '(' + mmbir + '){e<=' + str(mmbir_errors) + '}', ref)
  print(f'Found {len(output_list)} possible templates for {query}: {output_list}')

  if len(output_list) > 0:
    return(True)
  else:
    return(False)


def add_gene(row):

  data = EnsemblRelease(104)
  #print("loaded genome! #addgene")
  output = data.gene_names_at_locus(contig=row['chr'].strip("chr"), position=int(row["iBirStart"]))
  print(output)
  #output = str(output)
  out_list = []
  for gname in output:
    if gname != "":
      out_list.append(gname)
      #print("appended_gene!")

  return(out_list) #output

def add_exon(row):

  data = EnsemblRelease(104)
  #print("loaded genome! #addexon")
  output = data.exons_at_locus(contig=row['chr'].strip("chr"), position=int(row["iBirStart"]))
  #output = str(output)
  out_list = []
  for exon in output:
    if isinstance(exon, Exon):
      out_list.append(exon.exon_id)

  return(out_list) #output
            

def createAnnotatedOutput(f_name, path, output_f_name):
  chrList=createChrList(25)
  df = pd.DataFrame()

  for chr in chrList:
    try:
      with open(f"{path}{chr}/{f_name}") as f:
        lines = f.readlines()
    except:
      print(f"Couldn't read {path}{chr}/{f_name}. Make sure it's there.")
      lines = None
      continue

    if (type(lines) == None):
      print("Problem")

    else:
      print(f"hello {chr}")
      for i in range(len(lines)):
        if "###" in lines[i]:
          i-=1
          if "Microhomology Insertion:" in lines[i+11]:
            var0 = lines[i].strip()
            var1 = lines[i+1].strip()
            var2 = lines[i+2].strip("iBirStart:").strip()
            var3 = lines[i+3].strip("Consensus/cluster number:").strip()
            var4 = lines[i+4].strip("iDepth:").strip()
            var5 = lines[i+5].strip("sBir:").strip()
            var6 = lines[i+6].strip("sBirReversed:").strip()
            var7 = lines[i+7].strip()
            var8 = lines[i+8].strip("ref:").strip()
            var9 = lines[i+9].strip("bir:").strip()
            var10 = lines[i+10].strip("TEM:").strip()
            var11 = lines[i+11].strip("Microhomology Insertion:").strip()
            var12 = lines[i+12].strip("Microhomology template:").strip()
            var13 = lines[i+13].strip("readStart: ").strip()
            var14 = lines[i+14].strip("ref:").strip()
            var15 = lines[i+15].strip("bir:").strip()
            var16 = lines[i+16].strip()
            var17 = lines[i+17].strip()
            var18 = lines[i+18].strip()
            var19 = lines[i+19].strip()
          else:
            var0 = lines[i].strip()
            var1 = lines[i+1].strip()
            var2 = lines[i+2].strip("iBirStart:").strip()
            var3 = lines[i+3].strip("Consensus/cluster number:").strip()
            var4 = lines[i+4].strip("iDepth:").strip()
            var5 = lines[i+5].strip("sBir:").strip()
            var6 = lines[i+6].strip("sBirReversed:").strip()
            var7 = lines[i+7].strip()
            var8 = lines[i+8].strip("ref:").strip()
            var9 = lines[i+9].strip("bir:").strip()
            var10 = lines[i+10].strip("TEM:").strip()
            var11 = "Empty"
            var12 = "Empty"
            var13 = lines[i+11].strip("readStart: ").strip()
            var14 = lines[i+12].strip("ref:").strip()
            var15 = lines[i+13].strip("bir:").strip()
            var16 = lines[i+14].strip()
            var17 = lines[i+15].strip()
            var18 = lines[i+16].strip()

          event_dict = {
                        "chr":chr,
                        "iBirStart":var2,
                        "Consensus/cluster_number":var3,
                        "iDepth":var4,
                        "sBir":var5,
                        "sBirReversed":var6,
                        "ref":var8,
                        "bir":var9,
                        "TEM":var10,
                        "Microhomology_Insertion":var11,
                        "Microhomology_template":var12,
                        "readStart":var13,
                        "var16":var16,
                        "var17":var17,
                        "var18":var18
                        }
          event_df = pd.DataFrame.from_records([event_dict])
          df = pd.concat([df, event_df], ignore_index=True)

  print("Starting complexity")
  df[['ref_complexity_fail', 'ref_complexity_score', 'ref_complexity_tree']] = df.apply(lambda row: movingWindow(row.ref, complexity_threshold=0.2), axis=1, result_type ='expand')
  df[['bir_complexity_fail', 'bir_complexity_score', 'bir_complexity_tree']] = df.apply(lambda row: movingWindow(row.bir, complexity_threshold=0.2), axis=1, result_type ='expand')
  print("Starting homology check")
  df['homology_check_ref'] = df.apply(lambda row: verifyImperfectHomology(row.ref, row.sBir), axis=1)
  df['homology_check_bir'] = df.apply(lambda row: verifyImperfectHomology(row.bir, row.sBir), axis=1)
  print("Starting gene/exon annotation")
  df["genes"] = df.apply(lambda row: add_gene(row), axis=1)
  df["exones"] = df.apply(lambda row: add_exon(row), axis=1)

  df.to_csv(output_f_name, sep="\t",index=False)
  print(f"finished annotation... DF saved to {output_f_name}")


createAnnotatedOutput(f_name, path, sys.argv[2])
