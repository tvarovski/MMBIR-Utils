#!/usr/bin/python3

# this script is designed to be run from the command line and takes
# an annotated master_mmbir_table as input and annotates the sBir column
# with the complexity of the sBir

import os
import pandas as pd
from BioAid import findComplexity
# findComplexity is a function from the BioAid package and takes a sequence, a treeLevel, and a complexity_threshold as arguments
# and returns a list with the complexity, the treeLevel, and the sequence

def annotateSbirComplexity(df: pd.DataFrame, treeLevel: int, complexity_threshold: float) -> pd.DataFrame:
    '''
    Annotates the 'sBir' column of a DataFrame with a complexity measure.

    Args:
    - df: A pandas DataFrame that contains a column named 'sBir'.
    - treeLevel: An integer representing the level of the tree to consider when calculating complexity.
    - complexity_threshold: A float representing the threshold for complexity.

    Returns:
    - A pandas DataFrame with an additional column 'sBirComplexity' that contains the complexity measure for each 'sBir' value.

    This function applies the 'findComplexity' function to each value in the 'sBir' column of the DataFrame. The complexity measure is calculated based on the 'treeLevel' and 'complexity_threshold' parameters. The resulting complexity measure is stored in a new column 'sBirComplexity'.
    '''
    df["sBirComplexity"] = df["sBir"].apply(lambda x: findComplexity(x, treeLevel, complexity_threshold)[0])
    return df

def main(file_name: str, file_type: str = 'csv', tree_level: int = 2, complexity_threshold: float = 0.2):
    '''Main function to read a file, annotate the sBir column by complexity, and write the result to a new file.
    
    Args:
    - file_name: The name of the file to read.
    - file_type: The type of the file (default is 'csv'). Can be either 'csv' or 'tsv'.
    - tree_level: The level of the tree to consider when calculating complexity (default is 2).
    - complexity_threshold: The threshold for complexity (default is 0.2).
    
    Returns:
    - None
    '''

    print(f"{file_type} file option was selected.")
    sep = "\t" if file_type == "tsv" else ","
    df = pd.read_csv(file_name, sep=sep)

    df = annotateSbirComplexity(df, tree_level, complexity_threshold)

    base_name = os.path.splitext(file_name)[0]
    output_file = f"{base_name}_sBirAnn.{file_type}"

    df.to_csv(output_file, sep=sep, index=False)

if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2])
