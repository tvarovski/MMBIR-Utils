settings = {

    "TCGA-PROJECT": "KIRC", # e.g. "BRCA"
    "username": "twarowski", #e.g. jsmith123
    "consolidated_results_path": "consolidated_results_KIRC.tsv",
    "mmbir_fraction_high": 0.4,
    "mmbir_fraction_low": 0.4,
    "MMBIR_THRESHOLD_LOW": 150,
    "MMBIR_THRESHOLD_HIGH": 250,
    "min_concentration": 0,
    
    "ref_complexity_filter": True,
    "bir_complexity_filter": True,
    "ref_homology_check_filter": True,
    "bir_homology_check_filter": True,
    "exones_only": True,
    
    "cosmicdb_dir": 'tools/Cosmic.tsv',
    "outputs_path": "outputs",
    "outputs_path_raw": "outputs/raw",
    "outputs_path_filtered": "outputs/filtered",
    "expression_data_path_root": "/nfsscratch",
}