# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

import pandas as pd
import sys

def createChrList(chrNum):
  chrList = []
  for i in range(chrNum):
    chrList.append(f"chr{i+1}")
  return(chrList)

chrList=createChrList(25)


f_name = "final_bir_locs.txt"
#f_name = "non_commons_fromAll.txt"

path=sys.argv[1]
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

        df = df.append(event_dict, ignore_index=True)

df.to_csv(sys.argv[2], sep="\t",index=False)
