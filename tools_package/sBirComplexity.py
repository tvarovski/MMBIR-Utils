#!/usr/bin/python3

# this script is designed to be run from the command line and takes
# an annotated master_mmbir_table as input and annotates the sBir column
# with the complexity of the sBir

import sys
import pandas as pd
from BioAid import findComplexity
# findComplexity is a function from the BioAid package and takes a sequence, a treeLevel, and a complexity_threshold as arguments
# and returns a list with the complexity, the treeLevel, and the sequence

def annotateSbirComplexity(df, treeLevel, complexity_threshold):
    #annotate the sBir column by complexity, take only the first list element of output as the complexity
    df["sBirComplexity"] = df["sBir"].apply(lambda x: findComplexity(x, treeLevel, complexity_threshold)[0])
    return df

if __name__ == "__main__":
    # parameters

    treeLevel = 2
    complexity_threshold = 0.2

    # file can be either tsv or csv (default is csv)
    # specify with a command line argument
    file_name = sys.argv[1]
    file_type = sys.argv[2]

    # end parameters

    #read in the file
    if file_type == "tsv":
        print("tsv file option was selected.")
        df = pd.read_csv(file_name, sep="\t")
    else:
        print("csv file option was selected or no file type was specified, defaulting to csv.")
        df = pd.read_csv(file_name)


    #annotate the sBir column by complexity, take only the first list element of output as the complexity
    df = annotateSbirComplexity(df, treeLevel, complexity_threshold)

    #write to the file
    file_name = file_name.split(".")
    file_name = file_name[:-1]
    file_name = ".".join(file_name)

    if file_type == "tsv":
        df.to_csv(f"{file_name}_sBirAnn.tsv", sep="\t", index=False)
    else:
        df.to_csv(f"{file_name}_sBirAnn.csv", index=False)
