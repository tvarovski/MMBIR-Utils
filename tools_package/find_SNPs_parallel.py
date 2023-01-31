#this script takes a path to a text file and a number of processes as arguments
#and splits the file into chunks and runs the script on each chunk in parallel

import subprocess
import multiprocessing
import pandas as pd
import sys
from BioAid import split_df, join_slices, clean_up_slices


def runScriptSubprocess(df_slice):
    #this function runs the script on each chunk
    subprocess.call(['python', 'tools/find_SNPs.py', df_slice])


def runMainPool(results_df_path, num_processes=8):
    #this function splits the file into df split by chromosome (chr1, chr2, etc.)
    #and runs the script on each chromosome in parallel

    df = pd.read_csv(results_df_path)
    results_df_path_root = results_df_path[:-4]

    #split the df into n=num_processes slices and save each slice to a csv file
    print('Splitting df into slices...')
    split_df(df, num_processes, results_df_path_root)
    
    print(f'Running script in parallel in {num_processes} slices...')
    pool = multiprocessing.Pool(processes=num_processes)
    pool.map(runScriptSubprocess, [f"{results_df_path_root}_{slice}.csv" for slice in range(num_processes)])
    pool.close()
    pool.join()
    print(f'Finished {num_processes} processes...')

    #join the slices into one df
    print('Combining results from each slice...')
    join_slices(num_processes, results_df_path_root, f'{results_df_path_root}_SNP_anno_all.csv')

    #delete the slice files
    clean_up_slices(num_processes, results_df_path[:-4])


if __name__ == "__main__":

    try:
        results_df_path = sys.argv[1]
        num_processes = int(sys.argv[2])

        if results_df_path == "":
            print(f"no results_df_path given, using raw_mmbir_results_master_df.csv")
            results_df_path = "raw_mmbir_results_master_df.csv"
        
        if num_processes == "":
            print(f"no num_processes given, using 8")
            num_processes = 8

    except:
        print(f"no results_df_path or num_processes given, using raw_mmbir_results_master_df.csv and 8")
        exit()

    runMainPool(results_df_path, num_processes=num_processes)
