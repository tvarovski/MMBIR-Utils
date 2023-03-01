# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

import sys
import pandas as pd
import cancer_config as cfg
from tools import findCosmicGenes


if __name__ == "__main__":

    census_dir = cfg.settings['cosmicdb_dir']
    mmb_df_input_path = sys.argv[1]

    filter_dict={

        "ref_complexity_filter": cfg.settings["ref_complexity_filter"]
        "bir_complexity_filter": cfg.settings["bir_complexity_filter"]
        "ref_homology_check_filter": cfg.settings["ref_homology_check_filter"]
        "bir_homology_check_filter": cfg.settings["bir_homology_check_filter"]
        "exones_only": cfg.settings["exones_only"]

                  }

    findCosmicGenes(mmb_df_input_path, filter_dict)